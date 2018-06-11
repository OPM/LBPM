/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University

  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "common/StackTrace.h"

#include <algorithm>
#include <csignal>
#include <cstring>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>


#define perr std::cerr


// Detect the OS
// clang-format off
#if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 ) || defined( _MSC_VER )
    #define USE_WINDOWS
    #define NOMINMAX
#elif defined( __APPLE__ )
    #define USE_MAC
    #define USE_NM
#elif defined( __linux ) || defined( __linux__ ) || defined( __unix ) || defined( __posix )
    #define USE_LINUX
    #define USE_NM
#else
    #error Unknown OS
#endif
// clang-format on


// Include system dependent headers
// clang-format off
// Detect the OS and include system dependent headers
#ifdef USE_WINDOWS
    #include <windows.h>
    #include <dbghelp.h>
    #include <DbgHelp.h>
    #include <TlHelp32.h>
    #include <Psapi.h>
    #include <process.h>
    #include <stdio.h>
    #include <tchar.h>
    #pragma comment( lib, "version.lib" ) // for "VerQueryValue"
#else
    #include <dlfcn.h>
    #include <execinfo.h>
    #include <sched.h>
    #include <sys/time.h>
    #include <ctime>
    #include <unistd.h>
    #include <sys/syscall.h>
#endif
#ifdef USE_MAC
    #include <mach-o/dyld.h>
    #include <mach/mach.h>
    #include <sys/sysctl.h>
    #include <sys/types.h>
#endif
// clang-format on


#ifdef __GNUC__
#define USE_ABI
#include <cxxabi.h>
#endif


#ifndef NULL_USE
#define NULL_USE( variable )                 \
    do {                                     \
        if ( 0 ) {                           \
            char *temp = (char *) &variable; \
            temp++;                          \
        }                                    \
    } while ( 0 )
#endif


// Set the callstack signal
#ifdef SIGRTMIN
#define CALLSTACK_SIG SIGRTMIN + 4
#else
#define CALLSTACK_SIG SIGUSR1
#define SIGRTMIN SIGUSR1
#define SIGRTMAX SIGUSR1
#endif


// Helper thread
static std::shared_ptr<std::thread> globalMonitorThread;


// Utility to break a string by a newline
static inline std::vector<std::string> breakString( const std::string &str )
{
    std::vector<std::string> strvec;
    size_t i1 = 0;
    size_t i2 = std::min( str.find( '\n', i1 ), str.length() );
    while ( i1 < str.length() ) {
        strvec.push_back( str.substr( i1, i2 - i1 ) );
        i1 = i2 + 1;
        i2 = std::min( str.find( '\n', i1 ), str.length() );
    }
    return strvec;
}


// Function to replace all instances of a string with another
static inline void strrep( std::string &str, const std::string &s, const std::string &r )
{
    size_t i = 0;
    while ( i < str.length() ) {
        i = str.find( s, i );
        if ( i == std::string::npos ) {
            break;
        }
        str.replace( i, s.length(), r );
        i += r.length();
    }
}


// Utility to strip the path from a filename
static inline std::string stripPath( const std::string &filename )
{
    if ( filename.empty() )
        return std::string();
    int i = 0;
    for ( i = (int) filename.size() - 1; i >= 0 && filename[i] != 47 && filename[i] != 92; i-- ) {
    }
    i = std::max( 0, i + 1 );
    return filename.substr( i );
}


// Inline function to subtract two addresses returning the absolute difference
static inline void *subtractAddress( void *a, void *b )
{
    return reinterpret_cast<void *>(
        std::abs( reinterpret_cast<long long int>( a ) - reinterpret_cast<long long int>( b ) ) );
}


#ifdef USE_WINDOWS
static BOOL __stdcall readProcMem( HANDLE hProcess,
                                   DWORD64 qwBaseAddress,
                                   PVOID lpBuffer,
                                   DWORD nSize,
                                   LPDWORD lpNumberOfBytesRead )
{
    SIZE_T st;
    BOOL bRet = ReadProcessMemory( hProcess, (LPVOID) qwBaseAddress, lpBuffer, nSize, &st );
    *lpNumberOfBytesRead = (DWORD) st;
    return bRet;
}
static inline std::string getCurrentDirectory()
{
    char temp[1024] = { 0 };
    GetCurrentDirectoryA( sizeof( temp ), temp );
    return temp;
}
namespace StackTrace {
BOOL GetModuleListTH32( HANDLE hProcess, DWORD pid );
BOOL GetModuleListPSAPI( HANDLE hProcess );
DWORD LoadModule( HANDLE hProcess, LPCSTR img, LPCSTR mod, DWORD64 baseAddr, DWORD size );
void LoadModules();
}; // namespace StackTrace
#endif


// Functions to copy data
static inline char *copy_in( size_t N, const void *data, char *ptr )
{
    memcpy( ptr, data, N );
    return ptr + N;
}
static inline const char *copy_out( size_t N, void *data, const char *ptr )
{
    memcpy( data, ptr, N );
    return ptr + N;
}


/****************************************************************************
 *  Utility to call system command and return output                         *
 ****************************************************************************/
#ifdef USE_WINDOWS
#define popen _popen
#define pclose _pclose
#endif
std::string StackTrace::exec( const std::string &cmd, int &code )
{
    signal( SIGCHLD, SIG_DFL ); // Clear child exited
    FILE *pipe = popen( cmd.c_str(), "r" );
    if ( pipe == nullptr )
        return std::string();
    std::string result = "";
    result.reserve( 1024 );
    while ( !feof( pipe ) ) {
        char buffer[257];
        buffer[256] = 0;
        if ( fgets( buffer, 128, pipe ) != nullptr )
            result += buffer;
    }
    auto status = pclose( pipe );
    code        = WEXITSTATUS( status );
    return result;
}


/****************************************************************************
 *  stack_info                                                               *
 ****************************************************************************/
void StackTrace::stack_info::clear()
{
    address  = nullptr;
    address2 = nullptr;
    object.clear();
    function.clear();
    filename.clear();
    line = -1;
}
bool StackTrace::stack_info::operator==( const StackTrace::stack_info &rhs ) const
{
    if ( address == rhs.address )
        return true;
    if ( address2 == rhs.address2 && object == rhs.object )
        return true;
    return false;
}
bool StackTrace::stack_info::operator!=( const StackTrace::stack_info &rhs ) const
{
    return !operator==( rhs );
}
int StackTrace::stack_info::getAddressWidth() const
{
    auto addr = reinterpret_cast<unsigned long long int>( address );
    if ( addr <= 0xFFFF )
        return 4;
    if ( addr <= 0xFFFFFFFF )
        return 8;
    if ( addr <= 0xFFFFFFFFFFFF )
        return 12;
    return 16;
}
std::string
StackTrace::stack_info::print( int widthAddress, int widthObject, int widthFunction ) const
{
    char tmp1[64], tmp2[64];
    sprintf( tmp1, "0x%%0%illx:  ", widthAddress );
    sprintf( tmp2, tmp1, reinterpret_cast<unsigned long long int>( address ) );
    std::string stack( tmp2 );
    sprintf( tmp2, "%i", line );
    std::string line_str( tmp2 );
    size_t N = stack.length();
    stack += stripPath( object );
    stack.resize( std::max<size_t>( stack.size(), N + widthObject ), ' ' );
    N = stack.length() + 2;
    stack += "  " + function;
    if ( !filename.empty() && line > 0 ) {
        stack.resize( std::max<size_t>( stack.size(), N + widthFunction ), ' ' );
        stack += "  " + stripPath( filename ) + ":" + line_str;
    } else if ( !filename.empty() ) {
        stack.resize( std::max<size_t>( stack.size(), N + widthFunction ), ' ' );
        stack += "  " + stripPath( filename );
    } else if ( line > 0 ) {
        stack += " : " + line_str;
    }
    return stack;
}
size_t StackTrace::stack_info::size() const
{
    return 2 * sizeof( void * ) + 4 * sizeof( int ) + object.size() + function.size() +
           filename.size();
}
char *StackTrace::stack_info::pack( char *ptr ) const
{
    int Nobj  = object.size();
    int Nfun  = function.size();
    int Nfile = filename.size();
    ptr       = copy_in( sizeof( void * ), &address, ptr );
    ptr       = copy_in( sizeof( void * ), &address2, ptr );
    ptr       = copy_in( sizeof( int ), &Nobj, ptr );
    ptr       = copy_in( sizeof( int ), &Nfun, ptr );
    ptr       = copy_in( sizeof( int ), &Nfile, ptr );
    ptr       = copy_in( sizeof( int ), &line, ptr );
    ptr       = copy_in( Nobj, object.data(), ptr );
    ptr       = copy_in( Nfun, function.data(), ptr );
    ptr       = copy_in( Nfile, filename.data(), ptr );
    return ptr;
}
const char *StackTrace::stack_info::unpack( const char *ptr )
{
    int Nobj, Nfun, Nfile;
    ptr = copy_out( sizeof( void * ), &address, ptr );
    ptr = copy_out( sizeof( void * ), &address2, ptr );
    ptr = copy_out( sizeof( int ), &Nobj, ptr );
    ptr = copy_out( sizeof( int ), &Nfun, ptr );
    ptr = copy_out( sizeof( int ), &Nfile, ptr );
    ptr = copy_out( sizeof( int ), &line, ptr );
    object.resize( Nobj );
    function.resize( Nfun );
    filename.resize( Nfile );
    ptr = copy_out( Nobj, &object.front(), ptr );
    ptr = copy_out( Nfun, &function.front(), ptr );
    ptr = copy_out( Nfile, &filename.front(), ptr );
    return ptr;
}
std::vector<char> StackTrace::stack_info::packArray( const std::vector<stack_info> &data )
{
    size_t size = sizeof( int );
    for ( const auto &i : data )
        size += i.size();
    std::vector<char> vec( size, 0 );
    char *ptr = vec.data();
    int N     = data.size();
    ptr       = copy_in( sizeof( int ), &N, ptr );
    for ( const auto &i : data )
        ptr = i.pack( ptr );
    return vec;
}
std::vector<StackTrace::stack_info> StackTrace::stack_info::unpackArray( const char *ptr )
{
    int N;
    ptr = copy_out( sizeof( int ), &N, ptr );
    std::vector<stack_info> data( N );
    for ( auto &i : data )
        ptr = i.unpack( ptr );
    return data;
}
#ifdef USE_MPI
static std::vector<char> pack( const std::vector<std::vector<StackTrace::stack_info>> &data )
{
    size_t size = sizeof( int );
    for ( const auto &i : data ) {
        size += sizeof( int );
        for ( size_t j = 0; j < i.size(); j++ )
            size += i[j].size();
    }
    std::vector<char> out( size, 0 );
    char *ptr = out.data();
    int N     = data.size();
    ptr       = copy_in( sizeof( int ), &N, ptr );
    for ( int i = 0; i < N; i++ ) {
        int M = data[i].size();
        ptr   = copy_in( sizeof( int ), &M, ptr );
        for ( int j = 0; j < M; j++ )
            ptr = data[i][j].pack( ptr );
    }
    return out;
}
static std::vector<std::vector<StackTrace::stack_info>> unpack( const std::vector<char> &in )
{
    const char *ptr = in.data();
    int N;
    ptr = copy_out( sizeof( int ), &N, ptr );
    std::vector<std::vector<StackTrace::stack_info>> data( N );
    for ( int i = 0; i < N; i++ ) {
        int M;
        ptr = copy_out( sizeof( int ), &M, ptr );
        data[i].resize( M );
        for ( int j = 0; j < M; j++ )
            ptr = data[i][j].unpack( ptr );
    }
    return data;
}
#endif


/****************************************************************************
 *  multi_stack_info                                                         *
 ****************************************************************************/
StackTrace::multi_stack_info::multi_stack_info( const std::vector<stack_info> &rhs )
{
    operator=( rhs );
}
StackTrace::multi_stack_info &StackTrace::multi_stack_info::
operator=( const std::vector<stack_info> &rhs )
{
    clear();
    if ( rhs.empty() )
        return *this;
    N     = 1;
    stack = rhs[0];
    if ( rhs.size() > 1 )
        add( rhs.size() - 1, &rhs[1] );
    return *this;
}
void StackTrace::multi_stack_info::clear()
{
    N = 0;
    stack.clear();
    children.clear();
}
void StackTrace::multi_stack_info::print2( const std::string &prefix,
                                           int w[3],
                                           std::vector<std::string> &text ) const
{
    if ( stack == stack_info() ) {
        for ( const auto &child : children )
            child.print2( "", w, text );
        return;
    }
    std::string line = prefix + "[" + std::to_string( N ) + "] " + stack.print( w[0], w[1], w[2] );
    text.push_back( line );
    std::string prefix2 = prefix + "  ";
    for ( size_t i = 0; i < children.size(); i++ ) {
        const auto &child = children[i];
        std::vector<std::string> text2;
        child.print2( "", w, text2 );
        for ( size_t j = 0; j < text2.size(); j++ ) {
            std::string line = prefix2 + text2[j];
            if ( children.size() > 1 && j > 0 && i < children.size() - 1 )
                line[prefix2.size()] = '|';
            text.push_back( line );
        }
    }
}
std::vector<std::string> StackTrace::multi_stack_info::print( const std::string &prefix ) const
{
    std::vector<std::string> text;
    int w[3] = { 0 };
    w[0]     = getAddressWidth();
    w[1]     = getObjectWidth();
    w[2]     = getFunctionWidth();
    print2( prefix, w, text );
    return text;
}
int StackTrace::multi_stack_info::getAddressWidth() const
{
    int w = stack.getAddressWidth();
    for ( const auto &child : children )
        w = std::max( w, child.getAddressWidth() );
    return w;
}
int StackTrace::multi_stack_info::getObjectWidth() const
{
    int w = std::min<int>( stripPath( stack.object ).size() + 1, 20 );
    for ( const auto &child : children )
        w = std::max( w, child.getObjectWidth() );
    return w;
}
int StackTrace::multi_stack_info::getFunctionWidth() const
{
    int w = std::min<int>( stack.function.size() + 1, 40 );
    for ( const auto &child : children )
        w = std::max( w, child.getFunctionWidth() );
    return w;
}
void StackTrace::multi_stack_info::add( size_t len, const stack_info *stack )
{
    if ( len == 0 )
        return;
    const auto &s = stack[len - 1];
    for ( auto &i : children ) {
        if ( i.stack == s ) {
            i.N++;
            if ( len > 1 )
                i.add( len - 1, stack );
            return;
        }
    }
    children.resize( children.size() + 1 );
    children.back().N     = 1;
    children.back().stack = s;
    if ( len > 1 )
        children.back().add( len - 1, stack );
}


/****************************************************************************
 *  Function to find an entry                                                *
 ****************************************************************************/
template <class TYPE>
inline size_t findfirst( const std::vector<TYPE> &X, TYPE Y )
{
    if ( X.empty() )
        return 0;
    size_t lower = 0;
    size_t upper = X.size() - 1;
    if ( X[lower] >= Y )
        return lower;
    if ( X[upper] < Y )
        return upper;
    while ( ( upper - lower ) != 1 ) {
        size_t value = ( upper + lower ) / 2;
        if ( X[value] >= Y )
            upper = value;
        else
            lower = value;
    }
    return upper;
}


/****************************************************************************
 *  Function to get the executable name                                      *
 ****************************************************************************/
static char global_exe_name[1000] = { 0 };
static bool setGlobalExecutableName( char *exe )
{
    try {
#ifdef USE_LINUX
        auto *buf = new char[0x10000];
        int len   = ::readlink( "/proc/self/exe", buf, 0x10000 );
        if ( len != -1 ) {
            buf[len] = '\0';
            strcpy( exe, buf );
        }
        delete[] buf;
#elif defined( USE_MAC )
        uint32_t size = 0x10000;
        char *buf     = new char[size];
        memset( buf, 0, size );
        if ( _NSGetExecutablePath( buf, &size ) == 0 )
            strcpy( exe, buf );
        delete[] buf;
#elif defined( USE_WINDOWS )
        DWORD size = 0x10000;
        char *buf  = new char[size];
        memset( buf, 0, size );
        GetModuleFileName( nullptr, buf, size );
        strcpy( exe, buf );
        delete[] buf;
#endif
    } catch ( ... ) {
    }
    return true;
}
static bool global_exe_name_set = setGlobalExecutableName( global_exe_name );
std::string StackTrace::getExecutable()
{
    if ( !global_exe_name_set )
        global_exe_name_set = setGlobalExecutableName( global_exe_name );
    return std::string( global_exe_name );
}


/****************************************************************************
 * Function to get symbols for the executable from nm (if availible)         *
 * Note: this function maintains an internal cached copy to prevent          *
 *    exccessive calls to nm.  This function also uses a lock to ensure      *
 *    thread safety.                                                         *
 ****************************************************************************/
std::mutex getSymbols_mutex;
struct global_symbols_struct {
    std::vector<void *> address;
    std::vector<char> type;
    std::vector<std::string> obj;
    int error;
} global_symbols;
static const global_symbols_struct &getSymbols2()
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
#ifdef USE_LINUX
                sprintf( cmd, "nm -n --demangle %s", global_exe_name );
#elif defined( USE_MAC )
                sprintf( cmd, "nm -n %s | c++filt", global_exe_name );
#else
#error Unknown OS using nm
#endif
                int code;
                auto output = breakString( StackTrace::exec( cmd, code ) );
                for ( const auto &line : output ) {
                    if ( line.empty() )
                        continue;
                    if ( line[0] == ' ' )
                        continue;
                    auto *a = const_cast<char *>( line.c_str() );
                    char *b = strchr( a, ' ' );
                    if ( b == nullptr )
                        continue;
                    b[0] = 0;
                    b++;
                    char *c = strchr( b, ' ' );
                    if ( c == nullptr )
                        continue;
                    c[0] = 0;
                    c++;
                    char *d = strchr( c, '\n' );
                    if ( d )
                        d[0] = 0;
                    size_t add = strtoul( a, nullptr, 16 );
                    data.address.push_back( reinterpret_cast<void *>( add ) );
                    data.type.push_back( b[0] );
                    data.obj.emplace_back( c );
                }
            } catch ( ... ) {
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
int StackTrace::getSymbols( std::vector<void *> &address,
                            std::vector<char> &type,
                            std::vector<std::string> &obj )
{
    const global_symbols_struct &data = getSymbols2();
    address                           = data.address;
    type                              = data.type;
    obj                               = data.obj;
    return data.error;
}


/****************************************************************************
 *  Function to get call stack info                                          *
 ****************************************************************************/
#ifdef USE_MAC
static void *loadAddress( const std::string &object )
{
    static std::map<std::string, void *> obj_map;
    if ( obj_map.empty() ) {
        uint32_t numImages = _dyld_image_count();
        for ( uint32_t i = 0; i < numImages; i++ ) {
            const struct mach_header *header = _dyld_get_image_header( i );
            const char *name                 = _dyld_get_image_name( i );
            const char *p                    = strrchr( name, '/' );
            struct mach_header *address      = const_cast<struct mach_header *>( header );
            obj_map.insert( std::pair<std::string, void *>( p + 1, address ) );
            // printf("   module=%s, address=%p\n", p + 1, header);
        }
    }
    auto it       = obj_map.find( object );
    void *address = 0;
    if ( it != obj_map.end() ) {
        address = it->second;
    } else {
        it = obj_map.find( stripPath( object ) );
        if ( it != obj_map.end() )
            address = it->second;
    }
    // printf("%s: 0x%016llx\n",object.c_str(),address);
    return address;
}
static std::tuple<std::string, std::string, std::string, int> split_atos( const std::string &buf )
{
    if ( buf.empty() )
        return std::tuple<std::string, std::string, std::string, int>();
    // Get the function
    size_t index = buf.find( " (in " );
    if ( index == std::string::npos )
        return std::make_tuple(
            buf.substr( 0, buf.length() - 1 ), std::string(), std::string(), 0 );
    std::string fun = buf.substr( 0, index );
    std::string tmp = buf.substr( index + 5 );
    // Get the object
    index           = tmp.find( ')' );
    std::string obj = tmp.substr( 0, index );
    tmp             = tmp.substr( index + 1 );
    // Get the filename and line number
    size_t p1 = tmp.find( '(' );
    size_t p2 = tmp.find( ')' );
    tmp       = tmp.substr( p1 + 1, p2 - p1 - 1 );
    index     = tmp.find( ':' );
    std::string file;
    int line = 0;
    if ( index != std::string::npos ) {
        file = tmp.substr( 0, index );
        line = std::stoi( tmp.substr( index + 1 ) );
    } else if ( p1 != std::string::npos ) {
        file = tmp;
    }
    return std::make_tuple( fun, obj, file, line );
}
#endif
#ifdef USE_LINUX
using uint_p = uint64_t;
#elif defined( USE_MAC )
typedef unsigned long uint_p;
#endif
#if defined( USE_LINUX ) || defined( USE_MAC )
static inline std::string generateCmd( const std::string &s1,
                                       const std::string &s2,
                                       const std::string &s3,
                                       std::vector<void *> addresses,
                                       const std::string &s4 )
{
    std::string cmd = s1 + s2 + s3;
    for ( auto &addresse : addresses ) {
        char tmp[32];
        sprintf( tmp, "%lx ", reinterpret_cast<uint_p>( addresse ) );
        cmd += tmp;
    }
    cmd += s4;
    return cmd;
}
#endif
// clang-format off
static void getFileAndLineObject( std::vector<StackTrace::stack_info*> &info )
{
    if ( info.empty() )
        return;
    // This gets the file and line numbers for multiple stack lines in the same object
    #if defined( USE_LINUX )
        // Create the call command
        std::vector<void*> address_list(info.size(),nullptr);
        for (size_t i=0; i<info.size(); i++) {
            address_list[i] = info[i]->address;
            if ( info[i]->object.find( ".so" ) != std::string::npos )
                address_list[i] = info[i]->address2; 
            if ( info[i]->object.find( ".mexa64" ) != std::string::npos )
                address_list[i] = info[i]->address2; 
        }
        std::string cmd = generateCmd( "addr2line -C -e ", info[0]->object,
            " -f -i ", address_list, " 2> /dev/null" );
        // Get the function/line/file
        int code;
        auto cmd_output = StackTrace::exec( cmd, code );
        auto output = breakString( cmd_output );
        if ( output.size() != 2*info.size() )
            return;
        // Add the results to info
        for (size_t i=0; i<info.size(); i++) {
            // get function name
            if ( info[i]->function.empty() )
                info[i]->function = output[2*i+0];
            // get file and line
            const char *buf = output[2*i+1].c_str();
            if ( buf[0] != '?' && buf[0] != 0 ) {
                size_t j = 0;
                for ( j = 0; j < 4095 && buf[j] != ':'; j++ ) {
                }
                info[i]->filename = std::string( buf, j );
                info[i]->line     = atoi( &buf[j + 1] );
            }
        }
    #elif defined( USE_MAC )
        // Create the call command
        void* load_address = loadAddress( info[0]->object );
        if ( load_address == nullptr )
            return;
        std::vector<void*> address_list(info.size(),nullptr);
        for (size_t i=0; i<info.size(); i++)
            address_list[i] = info[i]->address;
        // Call atos to get the object info
        char tmp[64];
        sprintf( tmp, " -l %lx ", (uint_p) load_address );
        std::string cmd = generateCmd( "atos -o ", info[0]->object,
            tmp, address_list, " 2> /dev/null" );
        // Get the function/line/file
        int code;
        auto cmd_output = StackTrace::exec( cmd, code );
        auto output = breakString( cmd_output );
        if ( output.size() != info.size() )
            return;
        // Parse the output for function, file and line info
        for ( size_t i=0; i<info.size(); i++) {
            auto data = split_atos( output[i] );
            if ( info[i]->function.empty() )
                info[i]->function = std::get<0>(data);
            if ( info[i]->object.empty() )
                info[i]->object = std::get<1>(data);
            if ( info[i]->filename.empty() )
                info[i]->filename = std::get<2>(data);
            if ( info[i]->line==0 )
                info[i]->line = std::get<3>(data);
        }
    #endif
}
static void getFileAndLine( std::vector<StackTrace::stack_info> &info )
{
    // Build a list of stack elements for each object
    std::map<std::string,std::vector<StackTrace::stack_info*>> obj_map;
    for (auto & i : info) {
        auto& list = obj_map[i.object];
        list.emplace_back( &i );
    }
    // For each object, get the file/line numbers for all entries
    for ( auto& entry : obj_map ) 
        getFileAndLineObject( entry.second );
}
// Try to use the global symbols to decode info about the stack
static void getDataFromGlobalSymbols( StackTrace::stack_info &info )
{
    const global_symbols_struct &data = getSymbols2();
    if ( data.error == 0 ) {
        size_t index = findfirst( global_symbols.address, info.address );
        if ( index > 0 )
            info.object = global_symbols.obj[index - 1];
        else
            info.object = std::string(global_exe_name);
    }
}
static void signal_handler( int sig )
{
    printf("Signal caught acquiring stack (%i)\n",sig);
    StackTrace::setErrorHandlers( [](std::string,StackTrace::terminateType) { exit( -1 ); } );
}
StackTrace::stack_info StackTrace::getStackInfo( void *address )
{
    return getStackInfo( std::vector<void*>(1,address) )[0];
}
std::vector<StackTrace::stack_info> StackTrace::getStackInfo( const std::vector<void*>& address )
{
    // Temporarily handle signals to prevent recursion on the stack
    auto prev_handler = signal( SIGINT, signal_handler );
    // Get the detailed stack info
    std::vector<StackTrace::stack_info> info(address.size());
    try {
        #ifdef USE_WINDOWS
            IMAGEHLP_SYMBOL64 pSym[1024];
            memset( pSym, 0, sizeof( pSym ) );
            pSym->SizeOfStruct  = sizeof( IMAGEHLP_SYMBOL64 );
            pSym->MaxNameLength = 1024;

            IMAGEHLP_MODULE64 Module;
            memset( &Module, 0, sizeof( Module ) );
            Module.SizeOfStruct = sizeof( Module );

            HANDLE pid = GetCurrentProcess();

            for (size_t i=0; i<address.size(); i++) {
                info[i].address = address[i];
                DWORD64 address2 = reinterpret_cast<DWORD64>( address[i] );
                DWORD64 offsetFromSymbol;
                if ( SymGetSymFromAddr( pid, address2, &offsetFromSymbol, pSym ) != FALSE ) {
                    char name[8192]={0};
                    DWORD rtn = UnDecorateSymbolName( pSym->Name, name, sizeof(name)-1, UNDNAME_COMPLETE );
                    if ( rtn == 0 )
                        info[i].function = std::string(pSym->Name);
                    else
                        info[i].function = std::string(name);
                } else {
                    printf( "ERROR: SymGetSymFromAddr (%d,%p)\n", GetLastError(), address2 );
                }

                // Get line number
                IMAGEHLP_LINE64 Line;
                memset( &Line, 0, sizeof( Line ) );
                Line.SizeOfStruct = sizeof( Line );
                DWORD offsetFromLine;
                if ( SymGetLineFromAddr64( pid, address2, &offsetFromLine, &Line ) != FALSE ) {
                    info[i].line     = Line.LineNumber;
                    info[i].filename = std::string( Line.FileName );
                } else {
                    info[i].line     = 0;
                    info[i].filename = std::string();
                }

                // Get the object
                if ( SymGetModuleInfo64( pid, address2, &Module ) != FALSE ) {
                    //info[i].object = std::string( Module.ModuleName );
                    info[i].object = std::string( Module.LoadedImageName );
                    //info[i].baseOfImage = Module.BaseOfImage;
                }
            }
        #else
            for (size_t i=0; i<address.size(); i++) {
                info[i].address = address[i];
                #if defined(_GNU_SOURCE) || defined(USE_MAC)
                    Dl_info dlinfo;
                    if ( !dladdr( info[i].address, &dlinfo ) ) {
                        getDataFromGlobalSymbols( info[i] );
                        continue;
                    }
                    info[i].address2 = subtractAddress( info[i].address, dlinfo.dli_fbase );
                    info[i].object   = std::string( dlinfo.dli_fname );
                    #if defined( USE_ABI )
                        int status;
                        char *demangled = abi::__cxa_demangle( dlinfo.dli_sname, nullptr, nullptr, &status );
                        if ( status == 0 && demangled != nullptr ) {
                            info[i].function = std::string( demangled );
                        } else if ( dlinfo.dli_sname != nullptr ) {
                            info[i].function = std::string( dlinfo.dli_sname );
                        }
                        free( demangled );
                    #endif
                    if ( dlinfo.dli_sname != nullptr && info[i].function.empty() )
                        info[i].function = std::string( dlinfo.dli_sname );
                #else
                    getDataFromGlobalSymbols( info[i] );
                #endif
            }
            // Get the filename / line numbers for each item on the stack
            getFileAndLine( info );
        #endif
    } catch ( ... ) {
    }
    signal( SIGINT, prev_handler ) ;
    return info;
}


/****************************************************************************
*  Function to get the backtrace                                            *
****************************************************************************/
static int backtrace_thread( const std::thread::native_handle_type&, void**, size_t );
#if defined( USE_LINUX ) || defined( USE_MAC )
static int thread_backtrace_count;
static void* thread_backtrace[1000];
static std::mutex thread_backtrace_mutex;
static void _callstack_signal_handler( int, siginfo_t*, void* )
{
    thread_backtrace_count = backtrace_thread( StackTrace::thisThread(), thread_backtrace, 1000 );
}
#endif
static int backtrace_thread( const std::thread::native_handle_type& tid, void **buffer, size_t size )
{
    int count = 0;
    #if defined( USE_LINUX ) || defined( USE_MAC )
        // Get the trace
        if ( tid == pthread_self() ) {
            count = ::backtrace( buffer, size );
        } else {
            // Note: this will get the backtrace, but terminates the thread in the process!!!
            thread_backtrace_mutex.lock();
            struct sigaction sa;
            sigfillset(&sa.sa_mask);
            sa.sa_flags = SA_SIGINFO;
            sa.sa_sigaction = _callstack_signal_handler;
            sigaction(CALLSTACK_SIG, &sa, nullptr);
            thread_backtrace_count = -1;
            pthread_kill( tid, CALLSTACK_SIG );
            auto t1 = std::chrono::high_resolution_clock::now();
            auto t2 = std::chrono::high_resolution_clock::now();
            while ( thread_backtrace_count==-1 && std::chrono::duration<double>(t2-t1).count()<0.15 ) {
                std::this_thread::yield();
                t2 = std::chrono::high_resolution_clock::now();
            }
            count = std::max(thread_backtrace_count,0);
            memcpy( buffer, thread_backtrace, count*sizeof(void*) );
            thread_backtrace_count = -1;
            thread_backtrace_mutex.unlock();
        }
    #elif defined( USE_WINDOWS )
        #if defined(DBGHELP)

            // Load the modules for the stack trace
            LoadModules();

            // Initialize stackframe for first call
            ::CONTEXT context;
            memset( &context, 0, sizeof( context ) );
            context.ContextFlags = CONTEXT_FULL;
            RtlCaptureContext( &context );
            STACKFRAME64 frame; // in/out stackframe
            memset( &frame, 0, sizeof( frame ) );
            #ifdef _M_IX86
                DWORD imageType = IMAGE_FILE_MACHINE_I386;
                frame.AddrPC.Offset    = context.Eip;
                frame.AddrPC.Mode      = AddrModeFlat;
                frame.AddrFrame.Offset = context.Ebp;
                frame.AddrFrame.Mode   = AddrModeFlat;
                frame.AddrStack.Offset = context.Esp;
                frame.AddrStack.Mode   = AddrModeFlat;
            #elif _M_X64
                DWORD imageType = IMAGE_FILE_MACHINE_AMD64;
                frame.AddrPC.Offset    = context.Rip;
                frame.AddrPC.Mode      = AddrModeFlat;
                frame.AddrFrame.Offset = context.Rsp;
                frame.AddrFrame.Mode   = AddrModeFlat;
                frame.AddrStack.Offset = context.Rsp;
                frame.AddrStack.Mode   = AddrModeFlat;
            #elif _M_IA64
                DWORD imageType = IMAGE_FILE_MACHINE_IA64;
                frame.AddrPC.Offset     = context.StIIP;
                frame.AddrPC.Mode       = AddrModeFlat;
                frame.AddrFrame.Offset  = context.IntSp;
                frame.AddrFrame.Mode    = AddrModeFlat;
                frame.AddrBStore.Offset = context.RsBSP;
                frame.AddrBStore.Mode   = AddrModeFlat;
                frame.AddrStack.Offset  = context.IntSp;
                frame.AddrStack.Mode    = AddrModeFlat;
            #else
                #error "Platform not supported!"
            #endif

            auto pid = GetCurrentProcess();
            for ( int frameNum = 0; frameNum<1024; ++frameNum ) {
                BOOL rtn = StackWalk64( imageType, pid, tid, &frame, &context, readProcMem,
                                        SymFunctionTableAccess, SymGetModuleBase64, NULL );
                if ( !rtn ) {
                    printf( "ERROR: StackWalk64 (%p)\n", frame.AddrPC.Offset );
                    break;
                }
                if ( frame.AddrPC.Offset != 0 ) {
                    buffer[count] = reinterpret_cast<void*>( frame.AddrPC.Offset ) );
                    count++;
                }
                if ( frame.AddrReturn.Offset == 0 )
                    break;
            }
            SetLastError( ERROR_SUCCESS );
        #endif
    #else
        #warning Stack trace is not supported on this compiler/OS
    #endif
    return count;
}
std::vector<void*> StackTrace::backtrace( std::thread::native_handle_type tid )
{
    std::vector<void*> trace( 1000, nullptr );
    size_t count = backtrace_thread( tid, trace.data(), trace.size() );
    trace.resize(count);
    return trace;
}
std::vector<void*> StackTrace::backtrace()
{
    std::vector<void*> trace( 1000, nullptr );
    size_t count = backtrace_thread( thisThread(), trace.data(), trace.size() );
    trace.resize(count);
    return trace;
}
std::vector<std::vector<void *>> StackTrace::backtraceAll()
{
    // Get the list of threads
    auto threads = activeThreads( );
    // Get the backtrace of each thread
    std::vector<std::vector<void*>> trace(threads.size());
    size_t i = 0;
    for ( auto it=threads.begin(); i<threads.size(); i++, it++ ) {
        trace[i].resize(1000);
        size_t count = backtrace_thread( *it, trace[i].data(), trace[i].size() );
        trace[i].resize(count);
    }
    return trace;
}


/****************************************************************************
*  Function to get the list of all active threads                           *
****************************************************************************/
#if defined( USE_LINUX )
static std::thread::native_handle_type thread_handle;
static bool thread_id_finished;
static void _activeThreads_signal_handler( int )
{
    auto handle = StackTrace::thisThread( );
    thread_handle = handle;
    thread_id_finished = true;
}
static inline int get_tid( int pid, const std::string& line )
{
    char buf2[128]={0};
    int i1 = 0;
    while ( line[i1]==' ' && line[i1]!=0 ) { i1++; }
    int i2 = i1;
    while ( line[i2]!=' ' && line[i2]!=0 ) { i2++; }
    memcpy(buf2,&line[i1],i2-i1);
    buf2[i2-i1+1] = 0;
    int pid2 = atoi(buf2);
    if ( pid2 != pid )
        return -1;
    i1 = i2;
    while ( line[i1]==' ' && line[i1]!=0 ) { i1++; }
    i2 = i1;
    while ( line[i2]!=' ' && line[i2]!=0 ) { i2++; }
    memcpy(buf2,&line[i1],i2-i1);
    buf2[i2-i1+1] = 0;
    int tid = atoi(buf2);
    return tid;
}
#endif
std::thread::native_handle_type StackTrace::thisThread( )
{
    #if defined( USE_LINUX ) || defined( USE_MAC )
        return pthread_self();
    #elif defined( USE_WINDOWS )
        return GetCurrentThread();
    #else
        #warning Stack trace is not supported on this compiler/OS
        return std::thread::native_handle_type();
    #endif
}
std::set<std::thread::native_handle_type> StackTrace::activeThreads( )
{
    std::set<std::thread::native_handle_type> threads;
    #if defined( USE_LINUX )
        std::set<int> tid;
        int pid = getpid();
        char cmd[128];
        sprintf( cmd, "ps -T -p %i", pid );
        signal( SIGCHLD, SIG_DFL );     // Clear child exited
        int code;
        auto output = breakString( exec( cmd, code ) );
        for ( const auto& line : output ) {
            int tid2 = get_tid( pid, line );
            if ( tid2 != -1 )
                tid.insert( tid2 );
        }
        tid.erase( syscall(SYS_gettid) );
        signal( CALLSTACK_SIG, _activeThreads_signal_handler );
        for ( auto tid2 : tid ) {
            thread_backtrace_mutex.lock();
            thread_id_finished = false;
            thread_handle = thisThread();
            syscall( SYS_tgkill, pid, tid2, CALLSTACK_SIG );
            auto t1 = std::chrono::high_resolution_clock::now();
            auto t2 = std::chrono::high_resolution_clock::now();
            while ( !thread_id_finished && std::chrono::duration<double>(t2-t1).count()<0.1 ) {
                std::this_thread::yield();
                t2 = std::chrono::high_resolution_clock::now();
            }
            threads.insert( thread_handle );
            thread_backtrace_mutex.unlock();
        }
    #elif defined( USE_MAC )
        printf("activeThreads not finished\n");
    #elif defined( USE_WINDOWS )
        HANDLE hThreadSnap = CreateToolhelp32Snapshot( TH32CS_SNAPTHREAD, 0 ); 
        if( hThreadSnap != INVALID_HANDLE_VALUE ) {
            // Fill in the size of the structure before using it
            THREADENTRY32 te32
            te32.dwSize = sizeof(THREADENTRY32 );
            // Retrieve information about the first thread, and exit if unsuccessful
            if( !Thread32First( hThreadSnap, &te32 ) ) {
                printError( TEXT("Thread32First") );    // Show cause of failure
                CloseHandle( hThreadSnap );             // Must clean up the snapshot object!
                return( FALSE );
            }
            // Now walk the thread list of the system
            do { 
                if ( te32.th32OwnerProcessID == dwOwnerPID )
                    threads.insert( te32.th32ThreadID );
            } while( Thread32Next(hThreadSnap, &te32 ) );
            CloseHandle( hThreadSnap );                 // Must clean up the snapshot object!
        }
    #else
        #warning activeThreads is not yet supported on this compiler/OS
    #endif
    threads.insert( thisThread() );
    if ( globalMonitorThread )
        threads.erase( globalMonitorThread->native_handle() );
    return threads;
}
// clang-format on


/****************************************************************************
 *  Function to get the current call stack                                   *
 ****************************************************************************/
std::vector<StackTrace::stack_info> StackTrace::getCallStack()
{
    auto trace = StackTrace::backtrace();
    auto info  = getStackInfo( trace );
    return info;
}
std::vector<StackTrace::stack_info> StackTrace::getCallStack( std::thread::native_handle_type id )
{
    auto trace = StackTrace::backtrace( id );
    auto info  = getStackInfo( trace );
    return info;
}
static StackTrace::multi_stack_info
generateMultiStack( const std::vector<std::vector<void *>> &thread_backtrace )
{
    // Get the stack data for all pointers
    std::set<void *> addresses_set;
    for ( const auto &trace : thread_backtrace ) {
        for ( auto ptr : trace )
            addresses_set.insert( ptr );
    }
    std::vector<void *> addresses( addresses_set.begin(), addresses_set.end() );
    auto stack_data = StackTrace::getStackInfo( addresses );
    std::map<void *, StackTrace::stack_info> map_data;
    for ( size_t i = 0; i < addresses.size(); i++ )
        map_data.insert( std::make_pair( addresses[i], stack_data[i] ) );
    // Create the multi-stack trace
    StackTrace::multi_stack_info multistack;
    for ( const auto &trace : thread_backtrace ) {
        if ( trace.empty() )
            continue;
        // Create the stack for the given thread trace
        std::vector<StackTrace::stack_info> stack( trace.size() );
        for ( size_t i = 0; i < trace.size(); i++ )
            stack[i] = map_data[trace[i]];
        // Add the data to the multistack
        multistack.add( stack.size(), stack.data() );
    }
    return multistack;
}
StackTrace::multi_stack_info StackTrace::getAllCallStacks()
{
    // Get the backtrace of each thread
    auto thread_backtrace = backtraceAll();
    // Create the multi-stack strucutre
    auto stack = generateMultiStack( thread_backtrace );
    return stack;
}


/****************************************************************************
 *  Function to get system search paths                                      *
 ****************************************************************************/
std::string StackTrace::getSymPaths()
{
    std::string paths;
#ifdef USE_WINDOWS
    // Create the path list (seperated by ';' )
    paths = std::string( ".;" );
    paths.reserve( 1000 );
    // Add the current directory
    paths += getCurrentDirectory() + ";";
    // Now add the path for the main-module:
    char temp[1024];
    memset( temp, 0, sizeof( temp ) );
    if ( GetModuleFileNameA( nullptr, temp, sizeof( temp ) - 1 ) > 0 ) {
        for ( char *p = ( temp + strlen( temp ) - 1 ); p >= temp; --p ) {
            // locate the rightmost path separator
            if ( ( *p == '\\' ) || ( *p == '/' ) || ( *p == ':' ) ) {
                *p = 0;
                break;
            }
        }
        if ( strlen( temp ) > 0 ) {
            paths += temp;
            paths += ";";
        }
    }
    memset( temp, 0, sizeof( temp ) );
    if ( GetEnvironmentVariableA( "_NT_SYMBOL_PATH", temp, sizeof( temp ) - 1 ) > 0 ) {
        paths += temp;
        paths += ";";
    }
    memset( temp, 0, sizeof( temp ) );
    if ( GetEnvironmentVariableA( "_NT_ALTERNATE_SYMBOL_PATH", temp, sizeof( temp ) - 1 ) > 0 ) {
        paths += temp;
        paths += ";";
    }
    memset( temp, 0, sizeof( temp ) );
    if ( GetEnvironmentVariableA( "SYSTEMROOT", temp, sizeof( temp ) - 1 ) > 0 ) {
        paths += temp;
        paths += ";";
        // also add the "system32"-directory:
        paths += temp;
        paths += "\\system32;";
    }
    memset( temp, 0, sizeof( temp ) );
    if ( GetEnvironmentVariableA( "SYSTEMDRIVE", temp, sizeof( temp ) - 1 ) > 0 ) {
        paths += "SRV*;" + std::string( temp ) +
                 "\\websymbols*http://msdl.microsoft.com/download/symbols;";
    } else {
        paths += "SRV*c:\\websymbols*http://msdl.microsoft.com/download/symbols;";
    }
#endif
    return paths;
}


/****************************************************************************
 *  Load modules for windows                                                 *
 ****************************************************************************/
#ifdef USE_WINDOWS
BOOL StackTrace::GetModuleListTH32( HANDLE hProcess, DWORD pid )
{
    // CreateToolhelp32Snapshot()
    typedef HANDLE( __stdcall * tCT32S )( DWORD dwFlags, DWORD th32ProcessID );
    // Module32First()
    typedef BOOL( __stdcall * tM32F )( HANDLE hSnapshot, LPMODULEENTRY32 lpme );
    // Module32Next()
    typedef BOOL( __stdcall * tM32N )( HANDLE hSnapshot, LPMODULEENTRY32 lpme );

    // try both dlls...
    const TCHAR *dllname[] = { _T("kernel32.dll"), _T("tlhelp32.dll") };
    HINSTANCE hToolhelp    = nullptr;
    tCT32S pCT32S          = nullptr;
    tM32F pM32F            = nullptr;
    tM32N pM32N            = nullptr;

    HANDLE hSnap;
    MODULEENTRY32 me;
    me.dwSize = sizeof( me );

    for ( size_t i = 0; i < ( sizeof( dllname ) / sizeof( dllname[0] ) ); i++ ) {
        hToolhelp = LoadLibrary( dllname[i] );
        if ( hToolhelp == nullptr )
            continue;
        pCT32S = (tCT32S) GetProcAddress( hToolhelp, "CreateToolhelp32Snapshot" );
        pM32F  = (tM32F) GetProcAddress( hToolhelp, "Module32First" );
        pM32N  = (tM32N) GetProcAddress( hToolhelp, "Module32Next" );
        if ( ( pCT32S != nullptr ) && ( pM32F != nullptr ) && ( pM32N != nullptr ) )
            break; // found the functions!
        FreeLibrary( hToolhelp );
        hToolhelp = nullptr;
    }

    if ( hToolhelp == nullptr )
        return FALSE;

    hSnap = pCT32S( TH32CS_SNAPMODULE, pid );
    if ( hSnap == (HANDLE) -1 ) {
        FreeLibrary( hToolhelp );
        return FALSE;
    }

    bool keepGoing = !!pM32F( hSnap, &me );
    int cnt        = 0;
    while ( keepGoing ) {
        LoadModule( hProcess, me.szExePath, me.szModule, (DWORD64) me.modBaseAddr, me.modBaseSize );
        cnt++;
        keepGoing = !!pM32N( hSnap, &me );
    }
    CloseHandle( hSnap );
    FreeLibrary( hToolhelp );
    if ( cnt <= 0 )
        return FALSE;
    return TRUE;
}
DWORD StackTrace::LoadModule(
    HANDLE hProcess, LPCSTR img, LPCSTR mod, DWORD64 baseAddr, DWORD size )
{
    CHAR *szImg  = _strdup( img );
    CHAR *szMod  = _strdup( mod );
    DWORD result = ERROR_SUCCESS;
    if ( ( szImg == nullptr ) || ( szMod == nullptr ) ) {
        result = ERROR_NOT_ENOUGH_MEMORY;
    } else {
        if ( SymLoadModule( hProcess, 0, szImg, szMod, baseAddr, size ) == 0 )
            result = GetLastError();
    }
    ULONGLONG fileVersion = 0;
    if ( szImg != nullptr ) {
        // try to retrive the file-version:
        VS_FIXEDFILEINFO *fInfo = nullptr;
        DWORD dwHandle;
        DWORD dwSize = GetFileVersionInfoSizeA( szImg, &dwHandle );
        if ( dwSize > 0 ) {
            LPVOID vData = malloc( dwSize );
            if ( vData != nullptr ) {
                if ( GetFileVersionInfoA( szImg, dwHandle, dwSize, vData ) != 0 ) {
                    UINT len;
                    TCHAR szSubBlock[] = _T("\\");
                    if ( VerQueryValue( vData, szSubBlock, (LPVOID *) &fInfo, &len ) == 0 ) {
                        fInfo = nullptr;
                    } else {
                        fileVersion = ( (ULONGLONG) fInfo->dwFileVersionLS ) +
                                      ( (ULONGLONG) fInfo->dwFileVersionMS << 32 );
                    }
                }
                free( vData );
            }
        }

        // Retrive some additional-infos about the module
        IMAGEHLP_MODULE64 Module;
        Module.SizeOfStruct = sizeof( IMAGEHLP_MODULE64 );
        SymGetModuleInfo64( hProcess, baseAddr, &Module );
        LPCSTR pdbName = Module.LoadedImageName;
        if ( Module.LoadedPdbName[0] != 0 )
            pdbName = Module.LoadedPdbName;
    }
    if ( szImg != nullptr )
        free( szImg );
    if ( szMod != nullptr )
        free( szMod );
    return result;
}
BOOL StackTrace::GetModuleListPSAPI( HANDLE hProcess )
{
    DWORD cbNeeded;
    HMODULE hMods[1024];
    char tt[8192];
    char tt2[8192];
    if ( !EnumProcessModules( hProcess, hMods, sizeof( hMods ), &cbNeeded ) ) {
        return false;
    }
    if ( cbNeeded > sizeof( hMods ) ) {
        printf( "Insufficient memory allocated in GetModuleListPSAPI\n" );
        return false;
    }
    int cnt = 0;
    for ( DWORD i = 0; i < cbNeeded / sizeof( hMods[0] ); i++ ) {
        // base address, size
        MODULEINFO mi;
        GetModuleInformation( hProcess, hMods[i], &mi, sizeof( mi ) );
        // image file name
        tt[0] = 0;
        GetModuleFileNameExA( hProcess, hMods[i], tt, sizeof( tt ) );
        // module name
        tt2[0] = 0;
        GetModuleBaseNameA( hProcess, hMods[i], tt2, sizeof( tt2 ) );
        DWORD dwRes = LoadModule( hProcess, tt, tt2, (DWORD64) mi.lpBaseOfDll, mi.SizeOfImage );
        if ( dwRes != ERROR_SUCCESS )
            printf( "ERROR: LoadModule (%d)\n", dwRes );
        cnt++;
    }

    return cnt != 0;
}
void StackTrace::LoadModules()
{
    static bool modules_loaded = false;
    if ( !modules_loaded ) {
        modules_loaded = true;

        // Get the search paths for symbols
        std::string paths = StackTrace::getSymPaths();

        // Initialize the symbols
        if ( SymInitialize( GetCurrentProcess(), paths.c_str(), FALSE ) == FALSE )
            printf( "ERROR: SymInitialize (%d)\n", GetLastError() );

        DWORD symOptions = SymGetOptions();
        symOptions |= SYMOPT_LOAD_LINES | SYMOPT_FAIL_CRITICAL_ERRORS;
        symOptions     = SymSetOptions( symOptions );
        char buf[1024] = { 0 };
        if ( SymGetSearchPath( GetCurrentProcess(), buf, sizeof( buf ) ) == FALSE )
            printf( "ERROR: SymGetSearchPath (%d)\n", GetLastError() );

        // First try to load modules from toolhelp32
        BOOL loaded = StackTrace::GetModuleListTH32( GetCurrentProcess(), GetCurrentProcessId() );

        // Try to load from Psapi
        if ( !loaded )
            loaded = StackTrace::GetModuleListPSAPI( GetCurrentProcess() );
    }
}
#endif


/****************************************************************************
 *  Get the signal name                                                      *
 ****************************************************************************/
std::string StackTrace::signalName( int sig ) { return std::string( strsignal( sig ) ); }
std::vector<int> StackTrace::allSignalsToCatch()
{
    std::set<int> signals;
    for ( int i = 1; i < 32; i++ )
        signals.insert( i );
    for ( int i = SIGRTMIN; i <= SIGRTMAX; i++ )
        signals.insert( i );
    signals.erase( SIGKILL );
    signals.erase( SIGSTOP );
    return std::vector<int>( signals.begin(), signals.end() );
}
std::vector<int> StackTrace::defaultSignalsToCatch()
{
    auto tmp = allSignalsToCatch();
    std::set<int> signals( tmp.begin(), tmp.end() );
    signals.erase( SIGWINCH ); // Don't catch window changed by default
    signals.erase( SIGCONT );  // Don't catch continue by default
    return std::vector<int>( signals.begin(), signals.end() );
}


/****************************************************************************
 *  Set the signal handlers                                                  *
 ****************************************************************************/
static std::function<void( std::string, StackTrace::terminateType )> abort_fun;
static std::string rethrow()
{
    std::string last_message;
#ifdef USE_LINUX
    try {
        static int tried_throw = 0;
        if ( tried_throw == 0 ) {
            tried_throw = 1;
            throw;
        }
        // No active exception
    } catch ( const std::exception &err ) {
        // Caught a std::runtime_error
        last_message = err.what();
    } catch ( ... ) {
        // Caught an unknown exception
        last_message = "unknown exception occurred.";
    }
#endif
    return last_message;
}
static void term_func_abort( int sig )
{
    std::string msg( "Caught signal: " );
    msg += StackTrace::signalName( sig );
    abort_fun( msg, StackTrace::terminateType::signal );
}
static std::set<int> signals_set = std::set<int>();
static void term_func()
{
    std::string last_message = rethrow();
    StackTrace::clearSignals();
    abort_fun( "Unhandled exception:\n" + last_message, StackTrace::terminateType::exception );
}
void StackTrace::clearSignal( int sig )
{
    if ( signals_set.find( sig ) != signals_set.end() ) {
        signal( sig, SIG_DFL );
        signals_set.erase( sig );
    }
}
void StackTrace::clearSignals()
{
    for ( auto sig : signals_set )
        signal( sig, SIG_DFL );
    signals_set.clear();
}
void StackTrace::setSignals( const std::vector<int> &signals, void ( *handler )( int ) )
{
    for ( auto sig : signals ) {
        signal( sig, handler );
        signals_set.insert( sig );
    }
}
void StackTrace::setErrorHandlers(
    std::function<void( std::string, StackTrace::terminateType )> abort )
{
    abort_fun = abort;
    std::set_terminate( term_func );
    setSignals( defaultSignalsToCatch(), &term_func_abort );
    std::set_unexpected( term_func );
}


/****************************************************************************
 *  Global call stack functionallity                                         *
 ****************************************************************************/
#ifdef USE_MPI
static MPI_Comm globalCommForGlobalCommStack = MPI_COMM_NULL;
static bool stopGlobalMonitorThread          = false;
static void runGlobalMonitorThread()
{
    int rank = 0;
    int size = 1;
    MPI_Comm_size( globalCommForGlobalCommStack, &size );
    MPI_Comm_rank( globalCommForGlobalCommStack, &rank );
    while ( !stopGlobalMonitorThread ) {
        // Check for any messages
        int flag = 0;
        MPI_Status status;
        int err = MPI_Iprobe( MPI_ANY_SOURCE, 1, globalCommForGlobalCommStack, &flag, &status );
        if ( err != MPI_SUCCESS ) {
            printf( "Internal error in StackTrace::getGlobalCallStacks::runGlobalMonitorThread\n" );
            break;
        } else if ( flag != 0 ) {
            // We received a request
            int src_rank = status.MPI_SOURCE;
            int tag;
            MPI_Recv( &tag, 1, MPI_INT, src_rank, 1, globalCommForGlobalCommStack, &status );
            // Get a trace of all threads (except this)
            auto threads = StackTrace::activeThreads();
            threads.erase( StackTrace::thisThread() );
            if ( threads.empty() )
                continue;
            // Get the stack trace of each thread
            std::vector<std::vector<StackTrace::stack_info>> stack;
            for ( auto thread : threads )
                stack.push_back( StackTrace::getCallStack( thread ) );
            // Pack and send the data
            auto data = pack( stack );
            int count = data.size();
            MPI_Send( data.data(), count, MPI_CHAR, src_rank, tag, globalCommForGlobalCommStack );
        } else {
            // No requests recieved
            std::this_thread::sleep_for( std::chrono::milliseconds( 50 ) );
        }
    }
}
void StackTrace::globalCallStackInitialize( MPI_Comm comm )
{
#ifdef USE_MPI
    MPI_Comm_dup( comm, &globalCommForGlobalCommStack );
#endif
    stopGlobalMonitorThread = false;
    globalMonitorThread.reset( new std::thread( runGlobalMonitorThread ) );
}
void StackTrace::globalCallStackFinalize()
{
    stopGlobalMonitorThread = true;
    globalMonitorThread->join();
    globalMonitorThread.reset();
#ifdef USE_MPI
    if ( globalCommForGlobalCommStack != MPI_COMM_NULL )
        MPI_Comm_free( &globalCommForGlobalCommStack );
    globalCommForGlobalCommStack = MPI_COMM_NULL;
#endif
}
StackTrace::multi_stack_info StackTrace::getGlobalCallStacks()
{
    // Check if we properly initialized the comm
    if ( globalMonitorThread == nullptr ) {
        printf( "Warning: getGlobalCallStacks called without call to globalCallStackInitialize\n" );
        return getAllCallStacks();
    }
    if ( globalMonitorThread == nullptr ) {
        printf( "Warning: getGlobalCallStacks called without call to globalCallStackInitialize\n" );
        return getAllCallStacks();
    }
#ifdef USE_MPI
    int provided;
    MPI_Query_thread( &provided );
    if ( provided != MPI_THREAD_MULTIPLE ) {
        printf( "Warning: getGlobalCallStacks requires support for MPI_THREAD_MULTIPLE\n" );
        return getAllCallStacks();
    }
#endif
    if ( activeThreads().size() == 1 ) {
        printf( "Warning: getAllCallStacks not supported on this OS, defaulting to basic call "
                "stack\n" );
        return getAllCallStacks();
    }
    // Signal all processes that we want their stack for all threads
    int rank = 0;
    int size = 1;
    MPI_Comm_size( globalCommForGlobalCommStack, &size );
    MPI_Comm_rank( globalCommForGlobalCommStack, &rank );
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::uniform_int_distribution<> dis( 2, 0x7FFF );
    int tag = dis( gen );
    std::vector<MPI_Request> sendRequest( size );
    for ( int i = 0; i < size; i++ ) {
        if ( i == rank )
            continue;
        MPI_Isend( &tag, 1, MPI_INT, i, 1, globalCommForGlobalCommStack, &sendRequest[i] );
    }
    // Get the trace for the current process
    auto threads = StackTrace::activeThreads();
    StackTrace::multi_stack_info multistack;
    for ( auto thread : threads ) {
        auto stack = StackTrace::getCallStack( thread );
        multistack.add( stack.size(), stack.data() );
    }
    // Recieve the backtrace for all processes/threads
    int N_finished        = 1;
    auto start            = std::chrono::steady_clock::now();
    double time           = 0;
    const double max_time = 2.0 + size * 20e-3;
    while ( N_finished < size && time < max_time ) {
        int flag = 0;
        MPI_Status status;
        int err = MPI_Iprobe( MPI_ANY_SOURCE, tag, globalCommForGlobalCommStack, &flag, &status );
        if ( err != MPI_SUCCESS ) {
            printf( "Internal error in StackTrace::getGlobalCallStacks\n" );
            break;
        } else if ( flag != 0 ) {
            // We recieved a response
            int src_rank = status.MPI_SOURCE;
            int count;
            MPI_Get_count( &status, MPI_CHAR, &count );
            std::vector<char> data( count, 0 );
            MPI_Recv( data.data(),
                      count,
                      MPI_CHAR,
                      src_rank,
                      tag,
                      globalCommForGlobalCommStack,
                      &status );
            auto stack_list = unpack( data );
            for ( const auto &stack : stack_list )
                multistack.add( stack.size(), stack.data() );
            N_finished++;
        } else {
            auto stop = std::chrono::steady_clock::now();
            time      = std::chrono::duration_cast<std::chrono::seconds>( stop - start ).count();
            std::this_thread::yield();
        }
    }
    for ( int i = 0; i < size; i++ ) {
        if ( i == rank )
            continue;
        MPI_Request_free( &sendRequest[i] );
    }
    return multistack;
}
#else
void StackTrace::globalCallStackInitialize( MPI_Comm ) {}
void StackTrace::globalCallStackFinalize() {}
StackTrace::multi_stack_info StackTrace::getGlobalCallStacks() { return getAllCallStacks(); }
#endif


/****************************************************************************
 *  Cleanup the call stack                                                   *
 ****************************************************************************/
static inline size_t findMatching( const std::string &str, size_t pos )
{
    if ( str[pos] != '<' ) {
        perr << "Internal error string matching\n";
        perr << "   " << str << std::endl;
        perr << "   " << pos << std::endl;
        return pos;
    }
    size_t pos2 = pos + 1;
    int count   = 1;
    while ( count != 0 && pos2 < str.size() ) {
        if ( str[pos2] == '<' )
            count++;
        if ( str[pos2] == '>' )
            count--;
        pos2++;
    }
    return pos2;
}
void StackTrace::cleanupStackTrace( multi_stack_info &stack )
{
    auto it           = stack.children.begin();
    const size_t npos = std::string::npos;
    while ( it != stack.children.end() ) {
        auto &object      = it->stack.object;
        auto &function    = it->stack.function;
        auto &filename    = it->stack.filename;
        bool remove_entry = false;
        // Cleanup object and filename
        object   = stripPath( object );
        filename = stripPath( filename );
        // Remove callstack (and all children) for threads that are just contributing
        if ( function.find( "_callstack_signal_handler" ) != npos &&
             filename.find( "StackTrace.cpp" ) != npos ) {
            it = stack.children.erase( it );
            continue;
        }
        // Remove __libc_start_main
        if ( function.find( "__libc_start_main" ) != npos &&
             filename.find( "libc-start.c" ) != npos )
            remove_entry = true;
        // Remove backtrace_thread
        if ( function.find( "backtrace_thread" ) != npos &&
             filename.find( "StackTrace.cpp" ) != npos )
            remove_entry = true;
        // Remove __restore_rt
        if ( function.find( "__restore_rt" ) != npos && object.find( "libpthread" ) != npos )
            remove_entry = true;
        // Remove std::condition_variable::__wait_until_impl
        if ( function.find( "std::condition_variable::__wait_until_impl" ) != npos &&
             filename == "condition_variable" )
            remove_entry = true;
        // Remove std::_Function_handler<
        if ( function.find( "std::_Function_handler<" ) != npos && filename == "functional" )
            remove_entry = true;
        // Remove std::_Bind_simple<
        if ( function.find( "std::_Bind_simple<" ) != npos && filename == "functional" ) {
            auto pos     = function.find( "std::_Bind_simple<" );
            function     = function.substr( 0, pos ) + "std::_Bind_simple<...>(...)";
            remove_entry = true;
        }
        // Remove std::this_thread::__sleep_for
        if ( function.find( "std::this_thread::__sleep_for(" ) != npos &&
             object.find( "libstdc++" ) != npos )
            remove_entry = true;
        // Remove std::thread::_Impl
        if ( function.find( "std::thread::_Impl<" ) != npos && filename == "thread" )
            remove_entry = true;
        // Remove MATLAB internal routines
        if ( object == "libmwmcr.so" || object == "libmwm_lxe.so" || object == "libmwbridge.so" ||
             object == "libmwiqm.so" )
            remove_entry = true;
        // Remove the desired entry
        if ( remove_entry ) {
            if ( it->children.empty() ) {
                it = stack.children.erase( it );
                continue;
            } else if ( it->children.size() == 1 ) {
                *it = it->children[0];
                continue;
            }
        }
        // Cleanup template space
        strrep( function, " >", ">" );
        strrep( function, "< ", "<" );
        // Replace std::chrono::duration with abbriviated version
        if ( function.find( "std::chrono::duration<" ) != npos ) {
            strrep( function, "std::chrono::duration<long, std::ratio<1l, 1l> >", "ticks" );
            strrep( function,
                    "std::chrono::duration<long, std::ratio<1l, 1000000000l> >",
                    "nanoseconds" );
        }
        // Replace std::ratio with abbriviated version.
        if ( function.find( "std::ratio<" ) != npos ) {
            strrep( function, "std::ratio<1l, 1000000000000000000000000l>", "std::yocto" );
            strrep( function, "std::ratio<1l, 1000000000000000000000l>", "std::zepto" );
            strrep( function, "std::ratio<1l, 1000000000000000000l>", "std::atto" );
            strrep( function, "std::ratio<1l, 1000000000000000l>", "std::femto" );
            strrep( function, "std::ratio<1l, 1000000000000l>", "std::pico" );
            strrep( function, "std::ratio<1l, 1000000000l>", "std::nano" );
            strrep( function, "std::ratio<1l, 1000000l>", "std::micro" );
            strrep( function, "std::ratio<1l, 1000l>", "std::milli" );
            strrep( function, "std::ratio<1l, 100l>", "std::centi" );
            strrep( function, "std::ratio<1l, 10l>", "std::deci" );
            strrep( function, "std::ratio<1l, 1l>", "" );
            strrep( function, "std::ratio<10l, 1l>", "std::deca" );
            strrep( function, "std::ratio<60l, 1l>", "std::ratio<60>" );
            strrep( function, "std::ratio<100l, 1l>", "std::hecto" );
            strrep( function, "std::ratio<1000l, 1l>", "std::kilo" );
            strrep( function, "std::ratio<3600l, 1l>", "std::ratio<3600>" );
            strrep( function, "std::ratio<1000000l, 1l>", "std::mega" );
            strrep( function, "std::ratio<1000000000l, 1l>", "std::giga" );
            strrep( function, "std::ratio<1000000000000l, 1l>", "std::tera" );
            strrep( function, "std::ratio<1000000000000000l, 1l>", "std::peta" );
            strrep( function, "std::ratio<1000000000000000000l, 1l>", "std::exa" );
            strrep( function, "std::ratio<1000000000000000000000l, 1l>", "std::zetta" );
            strrep( function, "std::ratio<1000000000000000000000000l, 1l>", "std::yotta" );
            strrep( function, " >", ">" );
            strrep( function, "< ", "<" );
        }
        // Replace std::chrono::duration with abbriviated version.
        if ( function.find( "std::chrono::duration<" ) != npos ) {
            // clang-format off
            strrep( function, "std::chrono::duration<long, std::nano>", "std::chrono::nanoseconds" );
            strrep( function, "std::chrono::duration<long, std::micro>", "std::chrono::microseconds" );
            strrep( function, "std::chrono::duration<long, std::milli>", "std::chrono::milliseconds" );
            strrep( function, "std::chrono::duration<long>", "std::chrono::seconds" );
            strrep( function, "std::chrono::duration<long,>", "std::chrono::seconds" );
            strrep( function, "std::chrono::duration<long, std::ratio<60>>", "std::chrono::minutes" );
            strrep( function, "std::chrono::duration<long, std::ratio<3600>>", "std::chrono::hours" );
            strrep( function, " >", ">" );
            strrep( function, "< ", "<" );
            // clang-format on
        }
        // Replace std::this_thread::sleep_for with abbriviated version.
        if ( function.find( "::sleep_for<" ) != npos ) {
            strrep( function, "::sleep_for<long, std::nano>", "::sleep_for<nanoseconds>" );
            strrep( function, "::sleep_for<long, std::micro>", "::sleep_for<microseconds>" );
            strrep( function, "::sleep_for<long, std::milli>", "::sleep_for<milliseconds>" );
            strrep( function, "::sleep_for<long>", "::sleep_for<seconds>" );
            strrep( function, "::sleep_for<long,>", "::sleep_for<seconds>" );
            strrep( function, "::sleep_for<long, std::ratio<60>>", "::sleep_for<minutes>" );
            strrep( function, "::sleep_for<long, std::ratio<3600>>", "::sleep_for<hours>" );
            strrep( function,
                    "::sleep_for<nanoseconds>(std::chrono::nanoseconds",
                    "::sleep_for(std::chrono::nanoseconds" );
            strrep( function,
                    "::sleep_for<microseconds>(std::chrono::microseconds",
                    "::sleep_for(std::chrono::microseconds" );
            strrep( function,
                    "::sleep_for<milliseconds>(std::chrono::milliseconds",
                    "::sleep_for(std::chrono::milliseconds" );
            strrep( function,
                    "::sleep_for<seconds>(std::chrono::seconds",
                    "::sleep_for(std::chrono::seconds" );
            strrep( function,
                    "::sleep_for<milliseconds>(std::chrono::minutes",
                    "::sleep_for(std::chrono::milliseconds" );
            strrep( function,
                    "::sleep_for<milliseconds>(std::chrono::hours",
                    "::sleep_for(std::chrono::hours" );
        }
        // Replace std::basic_string with abbriviated version
        size_t pos = 0;
        while ( pos < function.size() ) {
            // Find next instance of std::basic_string
            const std::string match = "std::basic_string<";
            pos                     = function.find( match, pos );
            if ( pos == npos )
                break;
            // Find the matching >
            size_t pos1 = pos + match.size() - 1;
            size_t pos2 = findMatching( function, pos1 );
            if ( pos2 == pos1 )
                break;
            if ( function.substr( pos1 + 1, 4 ) == "char" )
                function.replace( pos, pos2 - pos, "std::string" );
            else if ( function.substr( pos1 + 1, 7 ) == "wchar_t" )
                function.replace( pos, pos2 - pos, "std::wstring" );
            else if ( function.substr( pos1 + 1, 8 ) == "char16_t" )
                function.replace( pos, pos2 - pos, "std::u16string" );
            else if ( function.substr( pos1 + 1, 8 ) == "char32_t" )
                function.replace( pos, pos2 - pos, "std::u32string" );
            pos++;
        }
        // Cleanup the children
        cleanupStackTrace( *it );
        ++it;
    }
}
