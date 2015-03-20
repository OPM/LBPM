# Create a shared_ptr.h file in the include directory that contains 
#    a shared_ptr class (hopefully typedef to a compiler basic)
# Arguements:
#    INSTALL_DIR - Directory to install shared_ptr.h
#    NAMESPACE - Namespace to contain the shared_ptr class (may be empty)
INCLUDE( CheckCXXSourceCompiles )
FUNCTION( CONFIGURE_SHARED_PTR INSTALL_DIR NAMESPACE )
    SET( CMAKE_REQUIRED_FLAGS ${CMAKE_CXX_FLAGS} )
    CHECK_CXX_SOURCE_COMPILES(
	    "   #include <memory>
            namespace ${NAMESPACE} { using std::shared_ptr; }
	        int main() {
	            ${NAMESPACE}::shared_ptr<int> ptr;
	            return 0;
	        }
	    "
	    MEMORY_SHARED_PTR )
    CHECK_CXX_SOURCE_COMPILES(
	    "   #include <memory>
            namespace ${NAMESPACE} { using std::tr1::shared_ptr; }
	        int main() {
	            ${NAMESPACE}::shared_ptr<int> ptr;
	            return 0;
	        }
	    "
	    MEMORY_TR1_SHARED_PTR )
    CHECK_CXX_SOURCE_COMPILES(
	    "   #include <tr1/memory>
            namespace  ${NAMESPACE} { using std::tr1::shared_ptr; }
	        int main() {
	            ${NAMESPACE}::shared_ptr<int> ptr;
	            return 0;
	        }
	    "
	    TR1_MEMORY_TR1_SHARED_PTR )
    GET_DIRECTORY_PROPERTY( dirs INCLUDE_DIRECTORIES )
    SET( CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}" )
    SET( CMAKE_REQUIRED_INCLUDES ${dirs} "${BOOST_INCLUDE}" )
    CHECK_CXX_SOURCE_COMPILES(
	    "   #include \"boost/shared_ptr.hpp\"
            namespace  ${NAMESPACE} { using boost::shared_ptr; }
	        int main() {
	            ${NAMESPACE}::shared_ptr<int> ptr;
	            return 0;
	        }
	    "
	    BOOST_SHARED_PTR )
    WRITE_DUMMY_SHARED_PTR( "${NAMESPACE}" "${CMAKE_CURRENT_BINARY_DIR}/tmp/dummy_shared_ptr.h" )
    CHECK_CXX_SOURCE_COMPILES(
	    "   #include <iostream>
	        #include \"${CMAKE_CURRENT_BINARY_DIR}/tmp/dummy_shared_ptr.h\"
	        int main() {
	            ${NAMESPACE}::shared_ptr<int> ptr;
	            return 0;
	        }
	    "
	    DUMMY_SHARED_PTR )
    IF ( NOT NAMESPACE )
        SET( NAMESPACE " " )
    ENDIF()
    IF ( BOOST_SHARED_PTR )
        FILE(WRITE  "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "#include \"boost/shared_ptr.hpp\"\n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "#include \"boost/weak_ptr.hpp\"\n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "#include \"boost/enable_shared_from_this.hpp\"\n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "namespace ${NAMESPACE} {\n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "    using boost::shared_ptr; \n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "    using boost::dynamic_pointer_cast; \n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "    using boost::const_pointer_cast; \n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "    using boost::weak_ptr; \n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "    using boost::enable_shared_from_this; \n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "}\n")
    ELSEIF ( MEMORY_SHARED_PTR )
        IF ( ${NAMESPACE} STREQUAL "std" )
            FILE(WRITE  "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "#include <memory>\n")
        ELSE()
            FILE(WRITE  "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "#include <memory>\n")
            FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "namespace ${NAMESPACE} {\n")
            FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "    using std::shared_ptr; \n")
            FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "    using std::dynamic_pointer_cast; \n")
            FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "    using std::const_pointer_cast; \n")
            FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "    using std::weak_ptr; \n")
            FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "    using std::enable_shared_from_this; \n")
            FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "}\n")
        ENDIF()
    ELSEIF ( MEMORY_TR1_SHARED_PTR )
        FILE(WRITE  "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "#include <memory>\n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "namespace ${NAMESPACE} {\n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "    using std::tr1::shared_ptr; \n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "    using std::tr1::dynamic_pointer_cast; \n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "    using std::tr1::const_pointer_cast; \n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "    using std::tr1::weak_ptr; \n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "    using std::tr1::enable_shared_from_this; \n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "}\n")
    ELSEIF ( TR1_MEMORY_TR1_SHARED_PTR )
        FILE(WRITE  "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "#include <tr1/memory>\n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "namespace ${NAMESPACE} {\n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "    using std::tr1::shared_ptr; \n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "    using std::tr1::dynamic_pointer_cast; \n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "    using std::tr1::const_pointer_cast; \n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "    using std::tr1::weak_ptr; \n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "    using std::tr1::enable_shared_from_this; \n")
        FILE(APPEND "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "}\n")
    ELSEIF ( DUMMY_SHARED_PTR ) 
        MESSAGE("Warning: No valid shared_ptr found, using dummy shared_ptr" )
        WRITE_DUMMY_SHARED_PTR( "${NAMESPACE}" "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" )
    ELSE()
        MESSAGE(FATAL_ERROR "No shared_ptr availible")
    ENDIF()
    EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_if_different 
        "${CMAKE_CURRENT_BINARY_DIR}/tmp/shared_ptr.h" "${INSTALL_DIR}/shared_ptr.h" )
ENDFUNCTION()


FUNCTION( WRITE_DUMMY_SHARED_PTR NAMESPACE FILENAME )
    FILE(WRITE  "${FILENAME}" "#ifndef DUMMY_SHARED_PTR_INC\n")
    FILE(APPEND "${FILENAME}" "#define DUMMY_SHARED_PTR_INC\n")
    FILE(APPEND "${FILENAME}" "namespace dummy {\n\n")
    FILE(APPEND "${FILENAME}" "template<class T> class shared_ptr {\n")
    FILE(APPEND "${FILENAME}" "public:\n")
    FILE(APPEND "${FILENAME}" "    shared_ptr( ): obj(NULL), count(NULL) {}\n")
    FILE(APPEND "${FILENAME}" "    shared_ptr( T *ptr ): obj(ptr), count(NULL) { if (ptr) { count = new int; (*count)=1; } } \n")
    FILE(APPEND "${FILENAME}" "    shared_ptr( const shared_ptr<T>& rhs ):  \n")
    FILE(APPEND "${FILENAME}" "        obj(rhs.get()), count(rhs.count) { if ( count!=NULL ) { ++(*count); } } \n")
    FILE(APPEND "${FILENAME}" "    template<class U> shared_ptr( const shared_ptr<U>& rhs ):  \n")
    FILE(APPEND "${FILENAME}" "        obj(rhs.get()), count(rhs.count) { if ( count!=NULL ) { ++(*count); } } \n")
    FILE(APPEND "${FILENAME}" "    shared_ptr& operator=( const shared_ptr<T>& rhs ) { obj=rhs.obj; count=rhs.count; ++(*count); return *this; } \n")
    FILE(APPEND "${FILENAME}" "    ~shared_ptr( ) { reset(); }\n")
    FILE(APPEND "${FILENAME}" "    void reset( T *ptr ) { reset(); obj=ptr; count=new int; (*count)=1; }\n")
    FILE(APPEND "${FILENAME}" "    void reset( void ) { \n")
    FILE(APPEND "${FILENAME}" "        if ( count!=NULL) { int tmp=--(*count); if ( tmp==0 ) { delete obj; delete count; } } \n")
    FILE(APPEND "${FILENAME}" "        obj=NULL; count=NULL; \n")
    FILE(APPEND "${FILENAME}" "    }\n")
    FILE(APPEND "${FILENAME}" "    T* get( ) const { return obj; } \n")
    FILE(APPEND "${FILENAME}" "    T* operator->( ) const { return obj; } \n")
    FILE(APPEND "${FILENAME}" "    const T& operator*( ) const { return *obj; } \n")
    FILE(APPEND "${FILENAME}" "    bool operator==( const T * rhs ) const { return obj==rhs; } \n")
    FILE(APPEND "${FILENAME}" "    bool operator!=( const T * rhs ) const { return obj!=rhs; } \n")
    FILE(APPEND "${FILENAME}" "protected:\n")
    FILE(APPEND "${FILENAME}" "    T *obj;\n")
    FILE(APPEND "${FILENAME}" "    volatile int *count;\n")
    FILE(APPEND "${FILENAME}" "template<class T1, class U> friend shared_ptr<T1> dynamic_pointer_cast( shared_ptr<U> const & );\n")
    FILE(APPEND "${FILENAME}" "template<class T1, class U> friend shared_ptr<T1> const_pointer_cast( shared_ptr<U> const & );\n")
    FILE(APPEND "${FILENAME}" "template<class Y> friend class shared_ptr;\n")
    FILE(APPEND "${FILENAME}" "};\n\n")
    FILE(APPEND "${FILENAME}" "template<class T, class U> shared_ptr<T> dynamic_pointer_cast( shared_ptr<U> const & rhs ) {\n")
    FILE(APPEND "${FILENAME}" "    T* obj = dynamic_cast<T*>(rhs.obj);\n")
    FILE(APPEND "${FILENAME}" "    shared_ptr<T> ptr;\n")
    FILE(APPEND "${FILENAME}" "    if ( obj!=NULL ) { ptr.obj = obj; ptr.count=rhs.count; ++(*ptr.count); }\n")
    FILE(APPEND "${FILENAME}" "    return ptr;\n}\n")
    FILE(APPEND "${FILENAME}" "template<class T, class U> shared_ptr<T> const_pointer_cast( shared_ptr<U> const & rhs ) {\n")
    FILE(APPEND "${FILENAME}" "    T* obj = const_cast<T*>(rhs.obj);\n")
    FILE(APPEND "${FILENAME}" "    shared_ptr<T> ptr;\n")
    FILE(APPEND "${FILENAME}" "    if ( obj!=NULL ) { ptr.obj = obj; ptr.count=rhs.count; ++(*ptr.count); }\n")
    FILE(APPEND "${FILENAME}" "    return ptr;\n}\n")
    FILE(APPEND "${FILENAME}" "\n} // namespace dummy\n")
    FILE(APPEND "${FILENAME}" "\n\n")
    FILE(APPEND "${FILENAME}" "namespace ${NAMESPACE} {\n")
    FILE(APPEND "${FILENAME}" "    using dummy::shared_ptr; \n")
    FILE(APPEND "${FILENAME}" "    using dummy::dynamic_pointer_cast; \n")
    FILE(APPEND "${FILENAME}" "    using dummy::const_pointer_cast; \n")
    FILE(APPEND "${FILENAME}" "}\n\n")
    FILE(APPEND "${FILENAME}" "#endif\n")
ENDFUNCTION()


