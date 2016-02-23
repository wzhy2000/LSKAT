// fm_new.h: interface for new/delete
//
//////////////////////////////////////////////////////////////////////

#if !defined(_FM_NEW_H_)
#define _FM_NEW_H_

#ifdef __cplusplus
extern "C" {
#endif

class CFmNewTemp {
public:
    void * allocate ( size_t size)
    {
		void* p = Calloc( size, char );
		return(p);
	}
    void deallocate ( void * p)
    {
		Free(p);
	}
} ;

void *operator new (std::size_t size, CFmNewTemp  & arena);
void operator delete (void * p, CFmNewTemp  & arena);


#ifdef __cplusplus
}
#endif

#endif


