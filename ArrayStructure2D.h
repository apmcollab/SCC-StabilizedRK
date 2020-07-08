//
//#####################################################################
//                   ArrayStructure2D.h  
//#####################################################################
//
// Class ArrayStructure2D
//
// Provides a "light weight" two dimensional array structure 
// with initialization capabilities and optional bounds checking. 
//
// This class can be used to provide an indexing structure 
// with bounds checking for an existing two dimensional array.
//
// The beginning index default is 0                    : (C convention)
// Data for the array is assumed to be stored by ROWS  : (C convention)
// Access using (*,*), e.g. A(i,j) for (i,j)th element.: (NOT C convention) 
//   
// If the class is created with a reference to an existing array, then
// deleting an instance of the class DOES NOT delete the data accessed
// by the instance.
//
// The copy constructor creates a duplicate instance. Deleting the copy 
// will not delete the data of the original. 
//
// The assignment operator creates a duplicate instance. Deleting the copy 
// will not delete the data of the original. 
//
// Specify the compiler directive _DEBUG when compiling to 
// force bounds checking. 
//
// This is a templated class, but the following defines are set at
// the end of the header to facilitate usage
//
// #define DoubleArrayStructure2D ArrayStructure2D<double>
// #define IntArrayStructure2D    ArrayStructure2D<int>
// #define FloatArrayStructure2D  ArrayStructure2D<float>
// #define LongArrayStructure2D   ArrayStructure2D<long>
// #define CharArrayStructure2D   ArrayStructure2D<char>
//
//#####################################################################
// Chris Anderson (C) UCLA                               Sept. 16, 1999
// Added void return for boundsCheck                     Jan. 24,2000
// Fixed memory leak for arrays initialized with 
// external data                                         June 25, 2000
//
// Added the assignment operator.                        Sept. 4, 2000
// Added initialize version of the copy constructor      May.  2, 2001
// Updated stream include to STL                         Jan.  23,2004
// Added initialize()                                     May   30,2005
//#####################################################################
//
#ifndef __ArrayStructure2D__
#define __ArrayStructure2D__

#ifdef  _DEBUG 
#include <iostream>
#include <cstdio>
using namespace std;
#endif

template <class T>
class ArrayStructure2D 
{
 
    public :

    ArrayStructure2D()
    {
    dataPtr       = 0;
    internalAlloc = 0;
    index1Size    = 0;
    index2Size    = 0;
    index1Begin   = 0;
    index2Begin   = 0;
    index1End     = 0;
    index2End     = 0;
    };


    ArrayStructure2D(long m, long n)
    {
    dataPtr       = 0;
    internalAlloc = 0; // will be reset in initialize
    initialize(m,n);
    };

    ArrayStructure2D(T* d, long m, long n)
    {
    internalAlloc = 0;
    initialize(d,m,n);
    };

    virtual ~ArrayStructure2D()
    {
    if(internalAlloc == 1) delete [] dataPtr;
    }

    ArrayStructure2D(const ArrayStructure2D& D)
    {
    if(D.dataPtr == 0)
    {
    dataPtr       = 0;
    internalAlloc = 0;
    index1Size    = 0;
    index2Size    = 0;
    index1Begin   = 0;
    index2Begin   = 0;
    index1End     = 0;
    index2End     = 0;
    return;
    }

    index1Size    = D.index1Size;
    index2Size    = D.index2Size;
    index1Begin   = D.index1Begin;
    index2Begin   = D.index2Begin;
    index1End     = D.index1End;
    index2End     = D.index2End;

	dataPtr       = new T[index1Size*index2Size];
    internalAlloc = 1;   
	long i;
	for(i = 0; i < index1Size*index2Size; i++) {dataPtr[i] = D.dataPtr[i];}
    };

    virtual void initialize()
    {
    if(internalAlloc == 1) delete [] dataPtr;
    dataPtr       = 0;
    internalAlloc = 0;
    index1Size    = 0;
    index2Size    = 0;
    index1Begin   = 0;
    index2Begin   = 0;
    index1End     = 0;
    index2End     = 0;
    };

    virtual void initialize(long m, long n)
    {
    if(internalAlloc == 1) 
    {
      if((index1Size != m)||(index2Size != n))
      {
        delete [] dataPtr;
        dataPtr = new T[m*n];
      }
    }
    else
    {
       if(dataPtr == 0)
       {
       dataPtr = new T[m*n];
       internalAlloc  = 1;
       }
    }
    index1Size    = m;
    index2Size    = n;
    index1Begin   = 0;
    index2Begin   = 0;
    index1End     = index1Begin + (index1Size - 1);
    index2End     = index2Begin + (index2Size - 1);
    };

    virtual void initialize(T* d, long m, long n)
    {
    if(internalAlloc == 1) delete [] dataPtr;
    dataPtr       = d;
    internalAlloc = 0;
    index1Size    = m;
    index2Size    = n;
    index1Begin   = 0;
    index2Begin   = 0;
    index1End     = index1Begin + (index1Size - 1);
    index2End     = index2Begin + (index2Size - 1);
    };

    virtual void initialize(const ArrayStructure2D& D)
    {

    if(internalAlloc == 1) 
    {
      if((index1Size != D.index1Size)||(index2Size != D.index2Size))
      {
        delete [] dataPtr;
        dataPtr = new T[(D.index1Size)*(D.index2Size)];
      }
    }
    else
    {
       if(dataPtr == 0)
       {
       dataPtr = new T[(D.index1Size)*(D.index2Size)];
       internalAlloc  = 1;
       }
    }

    index1Size    = D.index1Size;
    index2Size    = D.index2Size;
    index1Begin   = D.index1Begin;
    index2Begin   = D.index2Begin;
    index1End     = D.index1End;
    index2End     = D.index2End;

	long i;
	for(i = 0; i < index1Size*index2Size; i++) {dataPtr[i] = D.dataPtr[i];}
    };

    public :

#ifdef _DEBUG 
    T&  operator()(long i1, long i2)
    {
    boundsCheck(i1, index1Begin, index1End,1);
    boundsCheck(i2, index2Begin, index2End,2);
    return *(dataPtr +  (i2 - index2Begin) + (i1 - index1Begin)*index2Size);
    };

    const T&  operator()(long i1, long i2) const
    {
    boundsCheck(i1, index1Begin, index1End,1);
    boundsCheck(i2, index2Begin, index2End,2);
    return *(dataPtr +  (i2 - index2Begin) + (i1 - index1Begin)*index2Size);
    };
#else
    inline T&  operator()(long i1, long i2)
    {
    return *(dataPtr +  (i2 - index2Begin) + (i1 - index1Begin)*index2Size);
    };

    inline const T&  operator()(long i1, long i2) const
    {
    return *(dataPtr +  (i2 - index2Begin) + (i1 - index1Begin)*index2Size);
    };
#endif

void operator=(const ArrayStructure2D& D)
{

    #ifdef _DEBUG 
	if(index1Size != 0)
	{
    sizeCheck(this->index1Size,D.index1Size, this->index2Size,D.index2Size);
	}
	#endif

	if(index1Size*index2Size == 0)
	{
    initialize(D.index1Size,D.index2Size);
    
    index1Size    = D.index1Size;
    index2Size    = D.index2Size;
    index1Begin   = D.index1Begin;
    index2Begin   = D.index2Begin;
    index1End     = D.index1End;
    index2End     = D.index2End;
    }

    long i;
    for(i = 0; i < D.index1Size*D.index2Size; i++)
    {dataPtr[i] = D.dataPtr[i];}
    }

    public :

    T* getDataPointer(){return dataPtr;};

    const T* getDataPointer() const {return dataPtr;};

    void setIndex1Begin(long i) 
    {index1Begin = i; index1End   = index1Begin + (index1Size - 1);};

    void setIndex2Begin(long i)
    {index2Begin = i; index2End   = index2Begin + (index2Size - 1);};

    long getIndex1Begin() const {return index1Begin;}
    long getIndex2Begin() const {return index2Begin;}

    long getIndex1End() const {return index1End;}
    long getIndex2End() const {return index2End;}

    long getIndex1Size()  const {return index1Size;}
    long getIndex2Size()  const {return index2Size;}

 
	void setToValue(T d)
	{
    long i;
	for(i = 0; i < index1Size*index2Size; i++)
	{dataPtr[i] = d;}
	}
 
	void addValue(T d)
	{
    long i;
	for(i = 0; i < index1Size*index2Size; i++)
	{dataPtr[i] += d;}
	}
 
    protected :

    T*           dataPtr;     // data pointer
    long     index1Begin;     // coordinate 1 starting index
    long     index2Begin;     // coordinate 2 starting index
    long       index1End;     // coordinate 1 ending index
    long       index2End;     // coordinate 2 ending index
    long      index1Size;     // coordinate 1 size
    long      index2Size;     // coordinate 2 size
    int    internalAlloc;     // internal allocation flag 

#ifdef _DEBUG 
    static void boundsCheck(long i, long begin, long end, int coordinate)
    {
    if((i < begin)||(i  > end))
    {
    printf("Array index %d out of bounds \n",coordinate);
    printf("Offending index value %ld : Acceptable Range [%ld, %ld] \n",i, begin, end);
    }}
#else
static void boundsCheck(long, long, long, int){}
#endif

//
//###################################################################
// 
//###################################################################
//

#ifdef _DEBUG 
    static void sizeCheck(long Msize1, long Msize2, long Nsize1, long Nsize2)
    {
    if(Msize1 != Msize2)
    {
    printf("1st Dimension Sizes Are Incompatable  %ld != %ld \n" , Msize1, Msize2);
    }
	if(Nsize1 != Nsize2)
    {
    printf("2nd Dimension Sizes Are Incompatable  %ld != %ld \n" , Nsize1, Nsize2);
    }
    }
#else
static void sizeCheck(long, long, long, long){}
#endif
};


#define DoubleArrayStructure2D ArrayStructure2D<double>
#define IntArrayStructure2D    ArrayStructure2D<int>
#define FloatArrayStructure2D  ArrayStructure2D<float>
#define LongArrayStructure2D   ArrayStructure2D<long>
#define CharArrayStructure2D   ArrayStructure2D<char>
#endif 
