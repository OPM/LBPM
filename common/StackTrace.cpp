#include "StackTrace.h"

#include <iostream>
#include <sstream>
#include <cstring>
#include <algorithm>
#if __cplusplus > 199711L
    #include <mutex>
#endif


// Detect the OS and include system dependent headers
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64) || defined(_MSC_VER)
    // Note: windows has not been testeds
    #define USE_WINDOWS
    #include <iostream>
    #include <windows.h>
    #include <process.h>
    #include <stdio.h>   
    #include <tchar.h>
    #include <psapi.h>
    #include <DbgHelp.h>
    //#pragma comment(lib, psapi.lib) //added
    //#pragma comment(linker, /DEFAULTLIB:psapi.lib)
#elif defined(__APPLE__)
    #define USE_MAC
    #include <sys/time.h>
    #include <sys/sysctl.h>
    #include <signal.h>
    #include <execinfo.h>
    #include <dlfcn.h>
    #include <mach/mach.h>
    #include <sys/types.h>
    #include <sys/sysctl.h>
    #include <unistd.h>
    #include <sched.h>
#elif defined(__linux) || defined(__unix) || defined(__posix)
    #define USE_LINUX
    #define USE_NM
    #include <sys/time.h>
    #include <time.h>
    #include <execinfo.h>
    #include <dlfcn.h>
    #include <malloc.h>
    #include <unistd.h>
    #include <sched.h>
#else
    #error Unknown OS
#endif


#ifdef __GNUC__
    #define USE_ABI
    #include <cxxabi.h>
#endif

#ifndef NULL_USE
    #define NULL_USE(variable) do {                         \
        if(0) {char *temp = (char *)&variable; temp++;}     \
    }while(0)
#endif


// Utility to strip the path from a filename
inline std::string stripPath( const std::string& filename )
{
    if ( filename.empty() ) { return std::string(); }
    int i=0;
    for (i=(int)filename.size()-1; i>=0&&filename[i]!=47&&filename[i]!=92; i--) {}
    i = std::max(0,i+1);
    return filename.substr(i);
}


// Inline function to subtract two addresses returning the absolute difference
inline void* subtractAddress( void* a, void* b ) {
    return reinterpret_cast<void*>( std::abs(
        reinterpret_cast<long long int>(a)-reinterpret_cast<long long int>(b) ) );
}


/****************************************************************************
*  stack_info                                                               *
****************************************************************************/
std::string StackTrace::stack_info::print() const
{
    char tmp[32];
    sprintf(tmp,"0x%016llx:  ",reinterpret_cast<unsigned long long int>(address));
    std::string stack(tmp);
    sprintf(tmp,"%i",line);
    std::string line_str(tmp);
    stack += stripPath(object);
    stack.resize(std::max<size_t>(stack.size(),38),' ');
    stack += "  " + function;
    if ( !filename.empty() && line>0 ) {
        stack.resize(std::max<size_t>(stack.size(),70),' ');
        stack += "  " + stripPath(filename) + ":" + line_str;
    } else if ( !filename.empty() ) {
        stack.resize(std::max<size_t>(stack.size(),70),' ');
        stack += "  " + stripPath(filename);
    } else if ( line>0 ) {
        stack += " : " + line_str;
    }
    return stack;
}


/****************************************************************************
*  Function to find an entry                                                *
****************************************************************************/
template <class TYPE>
inline size_t findfirst( const std::vector<TYPE>& X, TYPE Y )
{
    if ( X.empty() )
        return 0;
    size_t lower = 0;
    size_t upper = X.size()-1;
    if ( X[lower] >= Y )
        return lower;
    if ( X[upper] < Y )
        return upper;
    while ( (upper-lower) != 1 ) {
        size_t value = (upper+lower)/2;
        if ( X[value] >= Y )
            upper = value;
        else
            lower = value;
    }
    return upper;
}


/****************************************************************************
* Function to get symbols for the executable from nm (if availible)         *
* Note: this function maintains an internal cached copy to prevent          *
*    exccessive calls to nm.  This function also uses a lock to ensure      *
*    thread safety.                                                         *
****************************************************************************/
#if __cplusplus <= 199711L
    class mutex_class {
      public:
        void lock() {}
        void unlock() {}
    };
    mutex_class getSymbols_mutex;
#else
    std::mutex getSymbols_mutex;
#endif
struct global_symbols_struct {
    std::vector<void*> address;
    std::vector<char> type;
    std::vector<std::string> obj;
    int error;
} global_symbols;
static std::string get_executable()
{
    std::string exe;
    try { 
        #ifdef USE_LINUX
            char *buf = new char[0x10000];
            int len = ::readlink("/proc/self/exe",buf,0x10000);
            if ( len!=-1 ) {
                buf[len] = '\0';
                exe = std::string(buf);
            }
        #endif
    } catch (...) {}
    return exe;
}
std::string global_exe_name = get_executable();
static const global_symbols_struct& getSymbols2(  )
{
    static bool loaded = false;
    static global_symbols_struct data;
    // Load the symbol tables if they have not been loaded
    if ( !loaded ) {
        getSymbols_mutex.lock();
        if ( !loaded ) {
            loaded = true;
            #ifdef USE_NM
                try { 
                    char cmd[1024];
                    sprintf(cmd,"nm --demangle --numeric-sort %s",global_exe_name.c_str());
                    FILE *in = popen(cmd,"r");
                    if ( in==NULL ) {
                        data.error = -2;
                        return data;
                    }
                    char *buf = new char[0x100000];
                    while ( fgets(buf,0xFFFFF,in)!=NULL ) {
                        if ( buf[0]==' ' || buf==NULL )
                            continue;
                        char *a = buf;
                        char *b = strchr(a,' ');  if (b==NULL) {continue;}  b[0] = 0;  b++;
                        char *c = strchr(b,' ');  if (c==NULL) {continue;}  c[0] = 0;  c++;
                        char *d = strchr(c,'\n');  if ( d ) { d[0]=0; }
                        size_t add = strtoul(a,NULL,16);
                        data.address.push_back( reinterpret_cast<void*>(add) );
                        data.type.push_back( b[0] );
                        data.obj.push_back( std::string(c) );
                    }
                    pclose(in);
                    delete [] buf;
                } catch (...) {
                    data.error = -3;
                }
                data.error = 0;
            #else
                data.error = -1;
            #endif
        }
        getSymbols_mutex.unlock();
    }
    return data;
}
int StackTrace::getSymbols( std::vector<void*>& address, std::vector<char>& type, 
    std::vector<std::string>& obj )
{
    const global_symbols_struct& data = getSymbols2();
    address = data.address;
    type = data.type;
    obj = data.obj;
    return data.error;
}


/****************************************************************************
*  Function to get the current call stack                                   *
****************************************************************************/
static void getFileAndLine( StackTrace::stack_info& info )
{
    #if defined(USE_LINUX) || defined(USE_MAC)
        void *address = info.address;
        if ( info.object.find(".so")!=std::string::npos )
            address = info.address2;
        char buf[4096];
        sprintf(buf, "addr2line -C -e %s -f -i %lx 2> /dev/null",
            info.object.c_str(),reinterpret_cast<unsigned long int>(address));
        FILE* f = popen(buf, "r");
        if (f == NULL)
            return;
        buf[4095] = 0;
        // get function name
        char *rtn = fgets(buf,4095,f);
        if ( info.function.empty() && rtn==buf ) {
            info.function = std::string(buf);
            info.function.resize(std::max<size_t>(info.function.size(),1)-1);
        }
        // get file and line
        rtn = fgets(buf,4095,f);
        if ( buf[0]!='?' && buf[0]!=0 && rtn==buf ) {
            size_t i = 0;
            for (i=0; i<4095 && buf[i]!=':'; i++) { }
            info.filename = std::string(buf,i);
            info.line = atoi(&buf[i+1]);
        }
        pclose(f);
    #endif
}

// Try to use the global symbols to decode info about the stack
static void getDataFromGlobalSymbols( StackTrace::stack_info& info )
{
    const global_symbols_struct& data = getSymbols2();
    if ( data.error==0 ) {
        size_t index = findfirst(global_symbols.address,info.address);
        if ( index > 0 )
            info.object = global_symbols.obj[index-1];
        else
            info.object = global_exe_name;
    }
}
StackTrace::stack_info StackTrace::getStackInfo( void* address )
{
    StackTrace::stack_info info;
    info.address = address;
    #ifdef _GNU_SOURCE
        Dl_info dlinfo;
        if ( !dladdr(address, &dlinfo) ) {
            getDataFromGlobalSymbols( info );
            getFileAndLine(info);
            return info;
        }
        info.address2 = subtractAddress(info.address,dlinfo.dli_fbase);
        info.object = std::string(dlinfo.dli_fname);
        #if defined(USE_ABI)
            int status;
            char *demangled = abi::__cxa_demangle(dlinfo.dli_sname,NULL,0,&status);
            if ( status == 0 && demangled!=NULL ) {
                info.function = std::string(demangled);
            } else if ( dlinfo.dli_sname!=NULL ) {
                info.function = std::string(dlinfo.dli_sname);
            }
            free(demangled);
        #else
            if ( dlinfo.dli_sname!=NULL )
                info.function = std::string(dlinfo.dli_sname);
        #endif
    #else
        getDataFromGlobalSymbols( info );
    #endif
    // Get the filename / line number
    getFileAndLine(info);
    return info;
}

std::vector<StackTrace::stack_info>  StackTrace::getCallStack()
{
    std::vector<StackTrace::stack_info>  stack_list;
    #if defined(USE_LINUX) || defined(USE_MAC)
        // Get the trace
        void *trace[100];
        memset(trace,0,100*sizeof(void*));
        int trace_size = backtrace(trace,100);
        stack_list.reserve(trace_size);
        for (int i=0; i<trace_size; ++i)
            stack_list.push_back(getStackInfo(trace[i]));
    #elif defined(USE_WINDOWS)
        #ifdef DBGHELP
            ::CONTEXT lContext;
            ::ZeroMemory( &lContext, sizeof( ::CONTEXT ) );
            ::RtlCaptureContext( &lContext );
            ::STACKFRAME64 lFrameStack;
            ::ZeroMemory( &lFrameStack, sizeof( ::STACKFRAME64 ) );
            lFrameStack.AddrPC.Offset = lContext.Rip;
            lFrameStack.AddrFrame.Offset = lContext.Rbp;
            lFrameStack.AddrStack.Offset = lContext.Rsp;
            lFrameStack.AddrPC.Mode = lFrameStack.AddrFrame.Mode = lFrameStack.AddrStack.Mode = AddrModeFlat;
            #ifdef _M_IX86
                DWORD MachineType = IMAGE_FILE_MACHINE_I386;
            #endif
            #ifdef _M_X64
                DWORD MachineType = IMAGE_FILE_MACHINE_AMD64;
            #endif
            #ifdef _M_IA64
                DWORD MachineType = IMAGE_FILE_MACHINE_IA64;
            #endif
            while ( 1 ) {
                int rtn = ::StackWalk64( MachineType, ::GetCurrentProcess(), ::GetCurrentThread(), 
                    &lFrameStack, MachineType == IMAGE_FILE_MACHINE_I386 ? 0 : &lContext,
                    NULL, &::SymFunctionTableAccess64, &::SymGetModuleBase64, NULL );
                if( !rtn )
                    break;
                if( lFrameStack.AddrPC.Offset == 0 )
                    break;
                ::MEMORY_BASIC_INFORMATION lInfoMemory;
                ::VirtualQuery( ( ::PVOID )lFrameStack.AddrPC.Offset, &lInfoMemory, sizeof( lInfoMemory ) );
                if ( lInfoMemory.Type==MEM_PRIVATE )
                    continue;
                ::DWORD64 lBaseAllocation = reinterpret_cast< ::DWORD64 >( lInfoMemory.AllocationBase );
                ::TCHAR lNameModule[ 1024 ];
                ::HMODULE hBaseAllocation = reinterpret_cast< ::HMODULE >( lBaseAllocation );
                ::GetModuleFileName( hBaseAllocation, lNameModule, 1024 );
                PIMAGE_DOS_HEADER lHeaderDOS = reinterpret_cast<PIMAGE_DOS_HEADER>( lBaseAllocation );
                if ( lHeaderDOS==NULL )
                    continue;
                PIMAGE_NT_HEADERS lHeaderNT = reinterpret_cast<PIMAGE_NT_HEADERS>( lBaseAllocation + lHeaderDOS->e_lfanew );
                PIMAGE_SECTION_HEADER lHeaderSection = IMAGE_FIRST_SECTION( lHeaderNT );
                ::DWORD64 lRVA = lFrameStack.AddrPC.Offset - lBaseAllocation;
                ::DWORD64 lNumberSection = ::DWORD64();
                ::DWORD64 lOffsetSection = ::DWORD64();
                for( int lCnt = ::DWORD64(); lCnt < lHeaderNT->FileHeader.NumberOfSections; lCnt++, lHeaderSection++ ) {
                    ::DWORD64 lSectionBase = lHeaderSection->VirtualAddress;
                    ::DWORD64 lSectionEnd = lSectionBase + std::max<::DWORD64>(
                        lHeaderSection->SizeOfRawData, lHeaderSection->Misc.VirtualSize );
                    if( ( lRVA >= lSectionBase ) && ( lRVA <= lSectionEnd ) ) {
                        lNumberSection = lCnt + 1;
                        lOffsetSection = lRVA - lSectionBase;
                        //break;
                    }
                }
                StackTrace::stack_info info;
                info.object = lNameModule;
                info.address = reinterpret_cast<void*>(lRVA);
                char tmp[20];
                sprintf(tmp,"0x%016llx",static_cast<unsigned long long int>(lOffsetSection));
                info.function = std::to_string(lNumberSection) + ":" + std::string(tmp);
                stack_list.push_back(info);
            }
        #endif
    #else
        #warning Stack trace is not supported on this compiler/OS
    #endif
    return stack_list;
}



