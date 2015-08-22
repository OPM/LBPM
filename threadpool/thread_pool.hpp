// This file contains the template functions for the thread pool
#ifndef included_ThreadPoolTmpl
#define included_ThreadPoolTmpl
#include "thread_pool.h"
#include <stdexcept>




/*! \addtogroup Macros
 *  @{
 */


/*! \def id = TPOOL_ADD_WORK(tpool,function,args,priority)
 *  \brief Add an item to the thread pool
 *  \details This a macro to automatically create and add a work item to 
 *      the thread pool.  
 *  \param tpool        Pointer to the thread pool to use
 *  \param function     Pointer to the function to use
 *  \param args         The arguments to pass to the function in the form (arg1,arg2,...)
 *  \param priority     Optional argument specifying the priority of the work item
 */
#define TPOOL_TUPLE_TO_SEQ(t) TPOOL_TUPLE_TO_SEQ_ ## II t
#define TPOOL_TUPLE_TO_SEQ_II(a,...) a,##__VA_ARGS__
#ifdef USE_WINDOWS
    #define TPOOL_GET_PRIORITY(a,N,c,...) N 
    #define TPOOL_ADD_WORK(TPOOL,FUNCTION,ARGS,...) \
        ThreadPool_add_work(TPOOL,FUNCTION,TPOOL_TUPLE_TO_SEQ(ARGS),TPOOL_GET_PRIORITY(0,__VA_ARGS__,0,0)+0)
#else
    #define TPOOL_GET_PRIORITY(_0,N,...) N 
    #define TPOOL_ADD_WORK(TPOOL,FUNCTION,ARGS,...) \
        ThreadPool_add_work(TPOOL,FUNCTION,TPOOL_TUPLE_TO_SEQ(ARGS),TPOOL_GET_PRIORITY(_0,##__VA_ARGS__,0))
#endif

/*! @} */

// \cond HIDDEN_SYMBOLS


// Specialization for no return argument
template <> 
class ThreadPool::WorkItemRet<void>: public ThreadPool::WorkItem {
public:
    virtual void run()=0;
    void get_results() { }
    virtual ~WorkItemRet() {}
};


// Final class for the work item
struct NULL_data {};
template < typename return_type, 
    typename  arg1=NULL_data,  typename  arg2=NULL_data,  typename  arg3=NULL_data,  typename  arg4=NULL_data, 
    typename  arg5=NULL_data,  typename  arg6=NULL_data,  typename  arg7=NULL_data,  typename  arg8=NULL_data,
    typename  arg9=NULL_data,  typename arg10=NULL_data,  typename arg11=NULL_data,  typename arg12=NULL_data, 
    typename arg13=NULL_data,  typename arg14=NULL_data,  typename arg15=NULL_data,  typename arg16=NULL_data,
    typename arg17=NULL_data,  typename arg18=NULL_data,  typename arg19=NULL_data,  typename arg20=NULL_data,
    typename arg21=NULL_data,  typename arg22=NULL_data,  typename arg23=NULL_data,  typename arg24=NULL_data > 
class WorkItemFull: public ThreadPool::WorkItemRet<return_type> 
{
private:
    int N;
    return_type (*routine)();
    arg1 x1;
    arg2 x2;
    arg3 x3;
    arg4 x4;
    arg5 x5;
    arg6 x6;
    arg7 x7;
    arg8 x8;
    arg9 x9;
    arg10 x10;
    arg11 x11;
    arg12 x12;
    arg13 x13;
    arg14 x14;
    arg15 x15;
    arg16 x16;
    arg17 x17;
    arg18 x18;
    arg19 x19;
    arg20 x20;
    arg21 x21;
    arg22 x22;
    arg23 x23;
    arg24 x24;
    WorkItemFull();
public:
    WorkItemFull( return_type (*routine2)() ):
        ThreadPool::WorkItemRet<return_type>(), N(0),
        routine(reinterpret_cast<return_type(*)()>(routine2)) { }
    WorkItemFull( return_type (*routine2)(arg1), arg1 y1):
        ThreadPool::WorkItemRet<return_type>(), N(1),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2), arg1 y1, arg2 y2 ):
        ThreadPool::WorkItemRet<return_type>(), N(2),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2,arg3), arg1 y1, arg2 y2, arg3 y3 ):
        ThreadPool::WorkItemRet<return_type>(), N(3),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2,arg3,arg4), arg1 y1, arg2 y2, arg3 y3, arg4 y4 ):
        ThreadPool::WorkItemRet<return_type>(), N(4),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2,arg3,arg4,arg5), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5 ):
        ThreadPool::WorkItemRet<return_type>(), N(5),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6 ):
        ThreadPool::WorkItemRet<return_type>(), N(6),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7 ):
        ThreadPool::WorkItemRet<return_type>(), N(7),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8 ):
        ThreadPool::WorkItemRet<return_type>(), N(8),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9 ):
        ThreadPool::WorkItemRet<return_type>(), N(9),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10 ):
        ThreadPool::WorkItemRet<return_type>(), N(10),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11 ):
        ThreadPool::WorkItemRet<return_type>(), N(11),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12 ):
        ThreadPool::WorkItemRet<return_type>(), N(12),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13 ):
        ThreadPool::WorkItemRet<return_type>(), N(13),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13, arg14 y14 ):
        ThreadPool::WorkItemRet<return_type>(), N(14),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13), x14(y14) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13, arg14 y14, arg15 y15 ):
        ThreadPool::WorkItemRet<return_type>(), N(15),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13), x14(y14), x15(y15) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13, arg14 y14, arg15 y15, arg16 y16 ):
        ThreadPool::WorkItemRet<return_type>(), N(16),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13), x14(y14), x15(y15), x16(y16) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13, arg14 y14, arg15 y15, arg16 y16, arg17 y17 ):
        ThreadPool::WorkItemRet<return_type>(), N(17),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13), x14(y14), x15(y15), x16(y16), x17(y17) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13, arg14 y14, arg15 y15, arg16 y16, arg17 y17, arg18 y18 ):
        ThreadPool::WorkItemRet<return_type>(), N(18),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13), x14(y14), x15(y15), x16(y16), x17(y17), x18(y18) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13, arg14 y14, arg15 y15, arg16 y16, arg17 y17, arg18 y18, arg19 y19 ):
        ThreadPool::WorkItemRet<return_type>(), N(19),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13), x14(y14), x15(y15), x16(y16), x17(y17), x18(y18), x19(y19) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13, arg14 y14, arg15 y15, arg16 y16, arg17 y17, arg18 y18, arg19 y19, arg20 y20 ):
        ThreadPool::WorkItemRet<return_type>(), N(20),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13), x14(y14), x15(y15), x16(y16), x17(y17), x18(y18), x19(y19), x20(y20) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13, arg14 y14, arg15 y15, arg16 y16, arg17 y17, arg18 y18, arg19 y19, arg20 y20, arg21 y21 ):
        ThreadPool::WorkItemRet<return_type>(), N(21),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13), x14(y14), x15(y15), x16(y16), x17(y17), x18(y18), x19(y19), x20(y20), x21(y21) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13, arg14 y14, arg15 y15, arg16 y16, arg17 y17, arg18 y18, arg19 y19, arg20 y20, arg21 y21, arg22 y22 ):
        ThreadPool::WorkItemRet<return_type>(), N(22),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13), x14(y14), x15(y15), x16(y16), x17(y17), x18(y18), x19(y19), x20(y20), x21(y21), x22(y22) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22,arg23), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13, arg14 y14, arg15 y15, arg16 y16, arg17 y17, arg18 y18, arg19 y19, arg20 y20, arg21 y21, arg22 y22, arg23 y23 ):
        ThreadPool::WorkItemRet<return_type>(), N(23),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13), x14(y14), x15(y15), x16(y16), x17(y17), x18(y18), x19(y19), x20(y20), x21(y21), x22(y22), x23(y23) { }
    WorkItemFull( return_type (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13, arg14 y14, arg15 y15, arg16 y16, arg17 y17, arg18 y18, arg19 y19, arg20 y20, arg21 y21, arg22 y22, arg23 y23, arg24 y24 ):
        ThreadPool::WorkItemRet<return_type>(), N(24),
        routine(reinterpret_cast<return_type(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13), x14(y14), x15(y15), x16(y16), x17(y17), x18(y18), x19(y19), x20(y20), x21(y21), x22(y22), x23(y23), x24(y24) { }
    void run() {
        ThreadPool::WorkItem::d_state = 1;
        if ( N==0 )
            this->d_result = reinterpret_cast<return_type(*)()>(routine)();
        else if ( N==1 )
            this->d_result = reinterpret_cast<return_type(*)(arg1)>(routine)(x1);
        else if ( N==2 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2)>(routine)(x1,x2);
        else if ( N==3 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2,arg3)>(routine)(x1,x2,x3);
        else if ( N==4 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2,arg3,arg4)>(routine)(x1,x2,x3,x4);
        else if ( N==5 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2,arg3,arg4,arg5)>(routine)(x1,x2,x3,x4,x5);
        else if ( N==6 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2,arg3,arg4,arg5,arg6)>(routine)(x1,x2,x3,x4,x5,x6);
        else if ( N==7 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7)>(routine)(x1,x2,x3,x4,x5,x6,x7);
        else if ( N==8 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8);
        else if ( N==9 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9);
        else if ( N==10 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10);
        else if ( N==11 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11);
        else if ( N==12 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12);
        else if ( N==13 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13);
        else if ( N==14 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14);
        else if ( N==15 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15);
        else if ( N==16 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16);
        else if ( N==17 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17);
        else if ( N==18 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18);
        else if ( N==19 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19);
        else if ( N==20 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20);
        else if ( N==21 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21);
        else if ( N==22 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22);
        else if ( N==23 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22,arg23)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23);
        else if ( N==24 )
            this->d_result = reinterpret_cast<return_type(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24);
        else
            throw std::exception(); // Internal error
        ThreadPool::WorkItem::d_state = 2;
    }
    virtual ~WorkItemFull() {}
};
template < 
typename  arg1,  typename  arg2,  typename  arg3,  typename  arg4,  typename  arg5,  
typename  arg6,  typename  arg7,  typename  arg8,  typename  arg9,  typename arg10,  
typename arg11,  typename arg12,  typename arg13,  typename arg14,  typename arg15,  
typename arg16,  typename arg17,  typename arg18,  typename arg19,  typename arg20,
typename arg21,  typename arg22,  typename arg23,  typename arg24 > 
class WorkItemFull<void,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24>: 
    public ThreadPool::WorkItemRet<void> 
{
private:
    int N;
    void (*routine)();
    arg1 x1;
    arg2 x2;
    arg3 x3;
    arg4 x4;
    arg5 x5;
    arg6 x6;
    arg7 x7;
    arg8 x8;
    arg9 x9;
    arg10 x10;
    arg11 x11;
    arg12 x12;
    arg13 x13;
    arg14 x14;
    arg15 x15;
    arg16 x16;
    arg17 x17;
    arg18 x18;
    arg19 x19;
    arg20 x20;
    arg21 x21;
    arg22 x22;
    arg23 x23;
    arg24 x24;
public:
    WorkItemFull( void (*routine2)() ):
        ThreadPool::WorkItemRet<void>(), N(0),
        routine(reinterpret_cast<void(*)()>(routine2)) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1), arg1 y1):
        ThreadPool::WorkItemRet<void>(), N(1),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2), arg1 y1, arg2 y2 ):
        ThreadPool::WorkItemRet<void>(), N(2),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2,arg3), arg1 y1, arg2 y2, arg3 y3 ):
        ThreadPool::WorkItemRet<void>(), N(3),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3){ ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2,arg3,arg4), arg1 y1, arg2 y2, arg3 y3, arg4 y4 ):
        ThreadPool::WorkItemRet<void>(), N(4),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2,arg3,arg4,arg5), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5 ):
        ThreadPool::WorkItemRet<void>(), N(5),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6 ):
        ThreadPool::WorkItemRet<void>(), N(6),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7 ):
        ThreadPool::WorkItemRet<void>(), N(7),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8 ):
        ThreadPool::WorkItemRet<void>(), N(8),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9 ):
        ThreadPool::WorkItemRet<void>(), N(9),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10 ):
        ThreadPool::WorkItemRet<void>(), N(10),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11 ):
        ThreadPool::WorkItemRet<void>(), N(11),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12 ):
        ThreadPool::WorkItemRet<void>(), N(12),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13 ):
        ThreadPool::WorkItemRet<void>(), N(13),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13, arg14 y14 ):
        ThreadPool::WorkItemRet<void>(), N(14),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13), x14(y14) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13, arg14 y14, arg15 y15 ):
        ThreadPool::WorkItemRet<void>(), N(15),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13), x14(y14), x15(y15) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13, arg14 y14, arg15 y15, arg16 y16 ):
        ThreadPool::WorkItemRet<void>(), N(16),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13), x14(y14), x15(y15), x16(y16) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13, arg14 y14, arg15 y15, arg16 y16, arg17 y17 ):
        ThreadPool::WorkItemRet<void>(), N(17),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13), x14(y14), x15(y15), x16(y16), x17(y17) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13, arg14 y14, arg15 y15, arg16 y16, arg17 y17, arg18 y18 ):
        ThreadPool::WorkItemRet<void>(), N(18),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13), x14(y14), x15(y15), x16(y16), x17(y17), x18(y18) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13, arg14 y14, arg15 y15, arg16 y16, arg17 y17, arg18 y18, arg19 y19 ):
        ThreadPool::WorkItemRet<void>(), N(19),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13), x14(y14), x15(y15), x16(y16), x17(y17), x18(y18), x19(y19) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13, arg14 y14, arg15 y15, arg16 y16, arg17 y17, arg18 y18, arg19 y19, arg20 y20 ):
        ThreadPool::WorkItemRet<void>(), N(20),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13), x14(y14), x15(y15), x16(y16), x17(y17), x18(y18), x19(y19), x20(y20) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13, arg14 y14, arg15 y15, arg16 y16, arg17 y17, arg18 y18, arg19 y19, arg20 y20, arg21 y21 ):
        ThreadPool::WorkItemRet<void>(), N(21),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13), x14(y14), x15(y15), x16(y16), x17(y17), x18(y18), x19(y19), x20(y20), x21(y21) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13, arg14 y14, arg15 y15, arg16 y16, arg17 y17, arg18 y18, arg19 y19, arg20 y20, arg21 y21, arg22 y22 ):
        ThreadPool::WorkItemRet<void>(), N(22),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13), x14(y14), x15(y15), x16(y16), x17(y17), x18(y18), x19(y19), x20(y20), x21(y21), x22(y22) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22,arg23), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13, arg14 y14, arg15 y15, arg16 y16, arg17 y17, arg18 y18, arg19 y19, arg20 y20, arg21 y21, arg22 y22, arg23 y23 ):
        ThreadPool::WorkItemRet<void>(), N(23),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13), x14(y14), x15(y15), x16(y16), x17(y17), x18(y18), x19(y19), x20(y20), x21(y21), x22(y22), x23(y23) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    WorkItemFull( void (*routine2)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24), arg1 y1, arg2 y2, arg3 y3, arg4 y4, arg5 y5, arg6 y6, arg7 y7, arg8 y8, arg9 y9, arg10 y10, arg11 y11, arg12 y12, arg13 y13, arg14 y14, arg15 y15, arg16 y16, arg17 y17, arg18 y18, arg19 y19, arg20 y20, arg21 y21, arg22 y22, arg23 y23, arg24 y24 ):
        ThreadPool::WorkItemRet<void>(), N(24),
        routine(reinterpret_cast<void(*)()>(routine2)),  x1(y1),  x2(y2),  x3(y3),  x4(y4),  x5(y5),  x6(y6),  x7(y7),  x8(y8),  x9(y9), x10(y10), x11(y11), x12(y12), x13(y13), x14(y14), x15(y15), x16(y16), x17(y17), x18(y18), x19(y19), x20(y20), x21(y21), x22(y22), x23(y23), x24(y24) {
        ThreadPool::WorkItem::d_state=0; ThreadPool::WorkItem::d_has_result=true; }
    void run() {
        ThreadPool::WorkItem::d_state = 1;
        if ( N==0 )
            reinterpret_cast<void(*)()>(routine)();
        else if ( N==1 )
            reinterpret_cast<void(*)(arg1)>(routine)(x1);
        else if ( N==2 )
            reinterpret_cast<void(*)(arg1,arg2)>(routine)(x1,x2);
        else if ( N==3 )
            reinterpret_cast<void(*)(arg1,arg2,arg3)>(routine)(x1,x2,x3);
        else if ( N==4 )
            reinterpret_cast<void(*)(arg1,arg2,arg3,arg4)>(routine)(x1,x2,x3,x4);
        else if ( N==5 )
            reinterpret_cast<void(*)(arg1,arg2,arg3,arg4,arg5)>(routine)(x1,x2,x3,x4,x5);
        else if ( N==6 )
            reinterpret_cast<void(*)(arg1,arg2,arg3,arg4,arg5,arg6)>(routine)(x1,x2,x3,x4,x5,x6);
        else if ( N==7 )
            reinterpret_cast<void(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7)>(routine)(x1,x2,x3,x4,x5,x6,x7);
        else if ( N==8 )
            reinterpret_cast<void(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8);
        else if ( N==9 )
            reinterpret_cast<void(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9);
        else if ( N==10 )
            reinterpret_cast<void(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10);
        else if ( N==11 )
            reinterpret_cast<void(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11);
        else if ( N==12 )
            reinterpret_cast<void(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12);
        else if ( N==13 )
            reinterpret_cast<void(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13);
        else if ( N==14 )
            reinterpret_cast<void(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14);
        else if ( N==15 )
            reinterpret_cast<void(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15);
        else if ( N==16 )
            reinterpret_cast<void(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16);
        else if ( N==17 )
            reinterpret_cast<void(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17);
        else if ( N==18 )
            reinterpret_cast<void(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18);
        else if ( N==19 )
            reinterpret_cast<void(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19);
        else if ( N==20 )
            reinterpret_cast<void(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20);
        else if ( N==21 )
            reinterpret_cast<void(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21);
        else if ( N==22 )
            reinterpret_cast<void(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22);
        else if ( N==23 )
            reinterpret_cast<void(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22,arg23)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23);
        else if ( N==24 )
            reinterpret_cast<void(*)(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24)>(routine)(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24);
        else
            throw std::exception(); // Internal error
        ThreadPool::WorkItem::d_state = 2;
    }
    virtual ~WorkItemFull() {}
};



// Function to get the returned function value
template <> 
inline bool ThreadPool::getFunctionRet<bool>( const ThreadPool::thread_id_t& id )  const
{    
    WorkItemRet<bool> *work = dynamic_cast<WorkItemRet<bool>*>(getFinishedWorkItem(id));
    bool rtn = false;
    if ( work != NULL )
        rtn = work->get_results();
    return rtn;
}
template <> 
inline char ThreadPool::getFunctionRet<char>( const ThreadPool::thread_id_t& id )  const
{    
    WorkItemRet<char> *work = dynamic_cast<WorkItemRet<char>*>(getFinishedWorkItem(id));
    char rtn = 0;
    if ( work != NULL )
        rtn = work->get_results();
    return rtn;
}
template <> 
inline int ThreadPool::getFunctionRet<int>( const ThreadPool::thread_id_t& id )  const
{    
    WorkItemRet<int> *work = dynamic_cast<WorkItemRet<int>*>(getFinishedWorkItem(id));
    int rtn = 0;
    if ( work != NULL )
        rtn = work->get_results();
    return rtn;
}
template <> 
inline float ThreadPool::getFunctionRet<float>( const ThreadPool::thread_id_t& id )  const
{    
    WorkItemRet<float> *work = dynamic_cast<WorkItemRet<float>*>(getFinishedWorkItem(id));
    float rtn = 0;
    if ( work != NULL )
        rtn = work->get_results();
    return rtn;
}
template <> 
inline double ThreadPool::getFunctionRet<double>( const ThreadPool::thread_id_t& id )  const
{    
    WorkItemRet<double> *work = dynamic_cast<WorkItemRet<double>*>(getFinishedWorkItem(id));
    double rtn = 0;
    if ( work != NULL )
        rtn = work->get_results();
    return rtn;
}
template <class return_type> 
inline return_type ThreadPool::getFunctionRet( const ThreadPool::thread_id_t& id )  const
{
    WorkItemRet<return_type> *work = dynamic_cast<WorkItemRet<return_type>*>(getFinishedWorkItem(id));
    return_type rtn;
    if ( work != NULL )
        rtn = work->get_results();
    return rtn;
}



// Functions create work items
template <class return_type> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)()) 
{
    WorkItemFull<return_type> *work;
    work = new WorkItemFull<return_type>(routine);
    return work;
}
template <class return_type> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(), void*) 
{
    WorkItemFull<return_type> *work;
    work = new WorkItemFull<return_type>(routine);
    return work;
}
template <class return_type, class type1> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1), type1 arg1) 
{
    WorkItemFull<return_type,type1> *work;
    work = new WorkItemFull<return_type,type1>(routine,arg1);
    return work;
}
template <class return_type, class type1, class type2> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2), type1 arg1, type2 arg2) 
{
    WorkItemFull<return_type,type1,type2> *work;
    work = new WorkItemFull<return_type,type1,type2>(routine,arg1,arg2);
    return work;
}
template <class return_type, class type1, class type2, class type3> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2,type3), 
    type1 arg1, type2 arg2, type3 arg3) 
{
    WorkItemFull<return_type,type1,type2,type3> *work;
    work = new WorkItemFull<return_type,type1,type2,type3>(routine,arg1,arg2,arg3);
    return work;
}
template <class return_type, class type1, class type2, class type3, class type4> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2,type3,type4), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4) 
{
    WorkItemFull<return_type,type1,type2,type3,type4> *work;
    work = new WorkItemFull<return_type,type1,type2,type3,type4>(routine,arg1,arg2,arg3,arg4);
    return work;
}
template <class return_type, class type1, class type2, class type3, class type4, class type5> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2,type3,type4,type5), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5) 
{
    WorkItemFull<return_type,type1,type2,type3,type4,type5> *work;
    work = new WorkItemFull<return_type,type1,type2,type3,type4,type5>(routine,arg1,arg2,arg3,arg4,arg5);
    return work;
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2,type3,type4,type5,type6), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6) 
{
    WorkItemFull<return_type,type1,type2,type3,type4,type5,type6> *work;
    work = new WorkItemFull<return_type,type1,type2,type3,type4,type5,type6>(routine,arg1,arg2,arg3,arg4,arg5,arg6);
    return work;
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2,type3,type4,type5,type6,type7), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7) 
{
    WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7> *work;
    work = new WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7>(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7);
    return work;
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8) 
{
    WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8> *work;
    work = new WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8>(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8);
    return work;
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9) 
{
    WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9> *work;
    work = new WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9>(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9);
    return work;
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10) 
{
    WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10> *work;
    work = new WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10>(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10);
    return work;
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11) 
{
    WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11> *work;
    work = new WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11>(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11);
    return work;
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12) 
{
    WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12> *work;
    work = new WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12>(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12);
    return work;
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13) 
{
    WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13> *work;
    work = new WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13>(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13);
    return work;
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13, class type14> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, type14 arg14) 
{
    WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14> *work;
    work = new WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14>(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14);
    return work;
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13, class type14, class type15> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, type14 arg14, type15 arg15) 
{
    WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15> *work;
    work = new WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15>(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15);
    return work;
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13, class type14, class type15, class type16> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, type14 arg14, type15 arg15, type16 arg16) 
{
    WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16> *work;
    work = new WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16>(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16);
    return work;
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13, class type14, class type15, class type16, class type17> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, type14 arg14, type15 arg15, type16 arg16, type17 arg17) 
{
    WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17> *work;
    work = new WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17>(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17);
    return work;
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13, class type14, class type15, class type16, class type17, class type18> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, type14 arg14, type15 arg15, type16 arg16, type17 arg17, type18 arg18) 
{
    WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18> *work;
    work = new WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18>(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18);
    return work;
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13, class type14, class type15, class type16, class type17, class type18, class type19> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, type14 arg14, type15 arg15, type16 arg16, type17 arg17, type18 arg18, type19 arg19) 
{
    WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19> *work;
    work = new WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19>(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19);
    return work;
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13, class type14, class type15, class type16, class type17, class type18, class type19, class type20> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19,type20), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, type14 arg14, type15 arg15, type16 arg16, type17 arg17, type18 arg18, type19 arg19, type20 arg20) 
{
    WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19,type20> *work;
    work = new WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19,type20>(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20);
    return work;
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13, class type14, class type15, class type16, class type17, class type18, class type19, class type20, class type21> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19,type20,type21), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, type14 arg14, type15 arg15, type16 arg16, type17 arg17, type18 arg18, type19 arg19, type20 arg20, type21 arg21) 
{
    WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19,type20,type21> *work;
    work = new WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19,type20,type21>(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21);
    return work;
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13, class type14, class type15, class type16, class type17, class type18, class type19, class type20, class type21, class type22> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19,type20,type21,type22), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, type14 arg14, type15 arg15, type16 arg16, type17 arg17, type18 arg18, type19 arg19, type20 arg20, type21 arg21, type22 arg22) 
{
    WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19,type20,type21,type22> *work;
    work = new WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19,type20,type21,type22>(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22);
    return work;
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13, class type14, class type15, class type16, class type17, class type18, class type19, class type20, class type21, class type22, class type23> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19,type20,type21,type22,type23), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, type14 arg14, type15 arg15, type16 arg16, type17 arg17, type18 arg18, type19 arg19, type20 arg20, type21 arg21, type22 arg22, type23 arg23) 
{
    WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19,type20,type21,type22,type23> *work;
    work = new WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19,type20,type21,type22,type23>(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22,arg23);
    return work;
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13, class type14, class type15, class type16, class type17, class type18, class type19, class type20, class type21, class type22, class type23, class type24> 
inline ThreadPool::WorkItem* ThreadPool_create_work(return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19,type20,type21,type22,type23,type24), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, type14 arg14, type15 arg15, type16 arg16, type17 arg17, type18 arg18, type19 arg19, type20 arg20, type21 arg21, type22 arg22, type23 arg23, type24 arg24) 
{
    WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19,type20,type21,type22,type23,type24> *work;
    work = new WorkItemFull<return_type,type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19,type20,type21,type22,type23,type24>(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24);
    return work;
}


// Functions to add work to the thread pool
template <class return_type> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(), int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work<return_type>(routine);
    return tpool->add_work( work, priority );
}
template <class return_type> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(), void*, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work<return_type>(routine);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1), type1 arg1, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work<return_type,type1>(routine,arg1);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2), type1 arg1, type2 arg2, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2, class type3> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2,type3), 
    type1 arg1, type2 arg2, type3 arg3, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2,arg3);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2, class type3, class type4> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2,type3,type4), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2,arg3,arg4);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2, class type3, class type4, class type5> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2,type3,type4,type5), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2,arg3,arg4,arg5);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2,type3,type4,type5,type6), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2,arg3,arg4,arg5,arg6);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2,type3,type4,type5,type6,type7), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13, class type14> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, type14 arg14, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13, class type14, class type15> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, type14 arg14, type15 arg15, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13, class type14, class type15, class type16> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, type14 arg14, type15 arg15, type16 arg16, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13, class type14, class type15, class type16, class type17> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, type14 arg14, type15 arg15, type16 arg16, type17 arg17, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13, class type14, class type15, class type16, class type17, class type18> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, type14 arg14, type15 arg15, type16 arg16, type17 arg17, type18 arg18, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13, class type14, class type15, class type16, class type17, class type18, class type19> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, type14 arg14, type15 arg15, type16 arg16, type17 arg17, type18 arg18, type19 arg19, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13, class type14, class type15, class type16, class type17, class type18, class type19, class type20> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19,type20), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, type14 arg14, type15 arg15, type16 arg16, type17 arg17, type18 arg18, type19 arg19, type20 arg20, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13, class type14, class type15, class type16, class type17, class type18, class type19, class type20, class type21> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19,type20,type21), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, type14 arg14, type15 arg15, type16 arg16, type17 arg17, type18 arg18, type19 arg19, type20 arg20, type21 arg21, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13, class type14, class type15, class type16, class type17, class type18, class type19, class type20, class type21, class type22> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19,type20,type21,type22), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, type14 arg14, type15 arg15, type16 arg16, type17 arg17, type18 arg18, type19 arg19, type20 arg20, type21 arg21, type22 arg22, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13, class type14, class type15, class type16, class type17, class type18, class type19, class type20, class type21, class type22, class type23> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19,type20,type21,type22,type23), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, type14 arg14, type15 arg15, type16 arg16, type17 arg17, type18 arg18, type19 arg19, type20 arg20, type21 arg21, type22 arg22, type23 arg23, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22,arg23);
    return tpool->add_work( work, priority );
}
template <class return_type, class type1, class type2, class type3, class type4, class type5, class type6, class type7, class type8, class type9, class type10, class type11, class type12, class type13, class type14, class type15, class type16, class type17, class type18, class type19, class type20, class type21, class type22, class type23, class type24> 
inline ThreadPool::thread_id_t ThreadPool_add_work(ThreadPool* tpool, return_type (*routine)(type1,type2,type3,type4,type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18,type19,type20,type21,type22,type23,type24), 
    type1 arg1, type2 arg2, type3 arg3, type4 arg4, type5 arg5, type6 arg6, type7 arg7, type8 arg8, type9 arg9, type10 arg10, type11 arg11, type12 arg12, type13 arg13, type14 arg14, type15 arg15, type16 arg16, type17 arg17, type18 arg18, type19 arg19, type20 arg20, type21 arg21, type22 arg22, type23 arg23, type24 arg24, int priority) 
{
    ThreadPool::WorkItem* work = ThreadPool_create_work(routine,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24);
    return tpool->add_work( work, priority );
}




/******************************************************************
* Inline function to return a unique id of the current thread     *
******************************************************************/
inline size_t ThreadPool::getThreadId()
{
    #ifdef USE_WINDOWS
        DWORD tmp_thread_id = GetCurrentThreadId();
        size_t thread_id = (size_t) tmp_thread_id;
    #elif defined(USE_LINUX) || defined(USE_MAC)
        pthread_t tmp_thread_id = pthread_self();
        size_t thread_id = (size_t) tmp_thread_id;
    #else 
        #error Not defined for this OS
    #endif
    return thread_id;
}


/******************************************************************
* Inline functions to wait for the work items to finish           *
******************************************************************/
inline int ThreadPool::wait( ThreadPool::thread_id_t id ) const
{
    bool finished;
    return wait_some( 1, &id, 1, &finished );
}
inline int ThreadPool::wait_any(size_t N_work, const ThreadPool::thread_id_t *ids) 
{
    bool* finished = new bool[N_work];
    int error = wait_some( N_work, ids, 1, finished );
    if ( error!=0 ) {
        delete [] finished;
        return error;
    }
    int index = -1;
    for (size_t i=0; i<N_work; i++) {
        if ( finished[i] ) {
            index = static_cast<int>(i);
            break;
        }
    }
    delete [] finished;
    return index;
}
inline int ThreadPool::wait_any(const std::vector<thread_id_t>& ids) const
{
    if ( ids.empty() )
        return 0;
    bool* finished = new bool[ids.size()];
    int error = wait_some( ids.size(), &ids[0], 1, finished );
    if ( error!=0 ) {
        delete [] finished;
        return error;
    }
    int index = -1;
    for (size_t i=0; i<ids.size(); i++) {
        if ( finished[i] ) {
            index = static_cast<int>(i);
            break;
        }
    }
    delete [] finished;
    return index;
}
inline int ThreadPool::wait_all(size_t N_work, const ThreadPool::thread_id_t *ids) const
{
    bool* finished = new bool[N_work];
    int error = wait_some( N_work, ids, N_work, finished );
    delete [] finished;
    return error;
}
inline int ThreadPool::wait_all(const std::vector<thread_id_t>& ids) const
{
    if ( ids.empty() )
        return 0;
    bool* finished = new bool[ids.size()];
    int error = wait_some( ids.size(), &ids[0], ids.size(), finished );
    delete [] finished;
    return error;
}


/******************************************************************
* Functions to add work items.                                    *
******************************************************************/
inline ThreadPool::thread_id_t ThreadPool::add_work( WorkItem *work, int priority ) 
{
    ThreadPool::thread_id_t id;
    add_work( 1, &work, &priority, &id );
    return id;
}
inline std::vector<ThreadPool::thread_id_t> ThreadPool::add_work( 
    const std::vector<ThreadPool::WorkItem*>& work, const std::vector<int>& priority ) 
{
    size_t N = work.size();
    if ( N==0 )
        return std::vector<ThreadPool::thread_id_t>();
    if ( priority.size()!=N && !priority.empty() )
        throw std::logic_error("size of work and priority do not match");
    const int* priority2 = NULL;
    if ( priority.empty() ) {
        priority2 = new int[N];
        memset( const_cast<int*>(priority2), 0, N*sizeof(int) );
    } else {
        priority2 = &priority[0];
    }
    std::vector<ThreadPool::thread_id_t> ids(N);
    add_work( N, const_cast<ThreadPool::WorkItem**>(&work[0]), priority2, &ids[0] );
    if ( priority.empty() )
        delete [] priority2;
    return ids;
}


/******************************************************************
* Class functions to for the thread id                            *
******************************************************************/
inline ThreadPool::thread_id_t::thread_id_t( ):
    d_id(0x0FFFFFFFFFFFFFFF), d_count(NULL), d_work(NULL)
{ }
inline ThreadPool::thread_id_t::~thread_id_t( )
{ 
    reset();
}
inline ThreadPool::thread_id_t::thread_id_t( const thread_id_t& rhs ):
    d_id(rhs.d_id), d_count(rhs.d_count), d_work(rhs.d_work)
{ 
    if ( d_count != NULL )
        AtomicOperations::atomic_increment(d_count);
}
inline ThreadPool::thread_id_t& ThreadPool::thread_id_t::operator=( const ThreadPool::thread_id_t& rhs ) volatile
{ 
    if (this == &rhs) // protect against invalid self-assignment
        return const_cast<ThreadPool::thread_id_t&>(*this);
    this->reset();
    d_id = rhs.d_id;
    d_count = rhs.d_count;
    d_work = rhs.d_work;
    if ( d_count != NULL )
        AtomicOperations::atomic_increment(d_count);
    return const_cast<ThreadPool::thread_id_t&>(*this);
}
#ifndef USE_WINDOWS
inline ThreadPool::thread_id_t::thread_id_t( const volatile ThreadPool::thread_id_t& rhs ):
    d_id(rhs.d_id), d_count(rhs.d_count), d_work(rhs.d_work)
{ 
    if ( d_count != NULL )
        AtomicOperations::atomic_increment(d_count);
}
inline ThreadPool::thread_id_t& ThreadPool::thread_id_t::operator=( const ThreadPool::thread_id_t& rhs )
{ 
    if (this == &rhs) // protect against invalid self-assignment
        return const_cast<ThreadPool::thread_id_t&>(*this);
    this->reset();
    d_id = rhs.d_id;
    d_count = rhs.d_count;
    d_work = rhs.d_work;
    if ( d_count != NULL )
        AtomicOperations::atomic_increment(d_count);
    return const_cast<ThreadPool::thread_id_t&>(*this);
}
inline ThreadPool::thread_id_t& ThreadPool::thread_id_t::operator=( const volatile ThreadPool::thread_id_t& rhs )
{ 
    if (this == &rhs) // protect against invalid self-assignment
        return const_cast<ThreadPool::thread_id_t&>(*this);
    this->reset();
    d_id = rhs.d_id;
    d_count = rhs.d_count;
    d_work = rhs.d_work;
    if ( d_count != NULL )
        AtomicOperations::atomic_increment(d_count);
    return const_cast<ThreadPool::thread_id_t&>(*this);
}
inline ThreadPool::thread_id_t& ThreadPool::thread_id_t::operator=( const volatile ThreadPool::thread_id_t& rhs ) volatile
{ 
    if (this == &rhs) // protect against invalid self-assignment
        return const_cast<ThreadPool::thread_id_t&>(*this);
    this->reset();
    d_id = rhs.d_id;
    d_count = rhs.d_count;
    d_work = rhs.d_work;
    if ( d_count != NULL )
        AtomicOperations::atomic_increment(d_count);
    return const_cast<ThreadPool::thread_id_t&>(*this);
}
#endif
inline void ThreadPool::thread_id_t::reset() volatile 
{ 
    if ( d_count != NULL ) {
        int count = AtomicOperations::atomic_decrement(d_count);
        if ( count == 0 ) {
            delete d_count;
            WorkItem* tmp = reinterpret_cast<ThreadPool::WorkItem*>(d_work);
            delete tmp;
        }
    }
    d_id = 0x0FFFFFFFFFFFFFFF;
    d_count = NULL;
    d_work = NULL;
}
inline void ThreadPool::thread_id_t::reset()
{ 
    if ( d_count != NULL ) {
        int count = AtomicOperations::atomic_decrement(d_count);
        if ( count == 0 ) {
            delete d_count;
            WorkItem* tmp = reinterpret_cast<ThreadPool::WorkItem*>(d_work);
            delete tmp;
        }
    }
    d_id = 0x0FFFFFFFFFFFFFFF;
    d_count = NULL;
    d_work = NULL;
}
inline void ThreadPool::thread_id_t::reset( int priority, size_t local_id, void* work )
{
    if ( d_count != NULL ) {
        int count = AtomicOperations::atomic_decrement(d_count);
        if ( count == 0 ) {
            delete d_count;
            WorkItem* tmp = reinterpret_cast<ThreadPool::WorkItem*>(d_work);
            delete tmp;
        }
    }
    // Create the id
    if ( sizeof(uint64)!=8 ) 
        throw std::logic_error("unsigned long long int must be 64 bits");        
    if ( priority<-127||priority>127 )
        throw std::logic_error("priority limited to +- 127");        
    if ( local_id > 0x00FFFFFFFFFFFFFF )
        throw std::logic_error("local id >= 2^56");
    char tmp1 = static_cast<char>(priority+128);
    unsigned char tmp2 = static_cast<unsigned char>(tmp1);
    if ( priority >= 0 )
        tmp2 |= 0x80;
    d_id = tmp2;
    d_id = (d_id<<56) + local_id;
    // Create the counter
    d_count = new AtomicOperations::int32_atomic;
    *d_count = 1;
    // Initialize the remaining data
    d_work = work;
}
inline size_t ThreadPool::thread_id_t::getLocalID(  ) const
{
    if ( d_id == 0x0FFFFFFFFFFFFFFF )
        return ~((size_t)0);
    unsigned long long int tmp = d_id&0x00FFFFFFFFFFFFFF;
    return static_cast<size_t>(tmp);
}
inline int ThreadPool::thread_id_t::getPriority(  ) const
{
    if ( d_id == 0x0FFFFFFFFFFFFFFF )
        return -128;
    unsigned long long int tmp = d_id>>56;
    return static_cast<int>(tmp)-128;
}
inline bool ThreadPool::thread_id_t::finished( ) const
{
    return d_id==0x0FFFFFFFFFFFFFFF ? true:reinterpret_cast<WorkItem*>(d_work)->d_state==2;
}


/******************************************************************
* This function checks if the id is valid                         *
******************************************************************/
#define MAXID32 0xFFFFFFFD
#define MAXID64 0x00FFFFFFFFFFFFFD
inline bool ThreadPool::isValid(const ThreadPool::thread_id_t& id) const
{
    size_t local_id = id.getLocalID();
    size_t next_id = d_id_assign-1;
    bool is_valid = true;
    if ( local_id==0 || !id.initialized() ) {
        // Invalid id
        is_valid = false;
    } else if ( d_id_assign==0 ) {
        // We ran out of thread ids
        throw std::logic_error("id space exhausted");
    } else if ( sizeof(size_t)==4 ) {
        // If we are using a 32-bit id, the work id is valid if it is <= 2^32-3 and > the next id
        if ( local_id>MAXID32 || local_id<=next_id )
            is_valid = false;
    } else if ( sizeof(size_t)==8 ) {
        // If we are using a 64-bit id, the work id is valid if it is <= 2^56-3, > the next id, and has an odd number of bits
        if ( local_id>MAXID64 || local_id<=next_id )
            is_valid = false;
    } else {
        // This is not a valid size
        throw std::logic_error("Error checking ids");
    }
    return is_valid;
}


// \endcond


#endif
