//
//#####################################################################
//                   ArrayStructure1D.h  
//#####################################################################
//
// Class ArrayStructure1D
//
// Provides a "light weight" one dimensional array structure 
// with initialization capabilities and optional bounds checking. 
//
// This class can be used to provide an indexing structure 
// with bounds checking for an existing one dimensional array. 
//
// The beginning index default is 0                : (C convention)
// Access using (*), e.g. A(i) for ith element.    : (Not C convention)
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
// This is a templated class, but the following defines are set at
// the end of the header to facilitate usage
//
// Specify the compiler directive _DEBUG when compiling to 
// force bounds checking. 
//
// #define DoubleArrayStructure1D ArrayStructure1D<double>
// #define IntArrayStructure1D    ArrayStructure1D<int>
// #define FloatArrayStructure1D  ArrayStructure1D<float>
// #define LongArrayStructure1D   ArrayStructure1D<long>
// #define CharArrayStructure1D   ArrayStructure1D<char>
//
//#####################################################################
// Chris Anderson (C) UCLA                               Sept 16, 1999
// Added void return for boundsCheck                     Jan. 24, 2000
// Fixed memory leak for arrays initialized with 
// external data                                         June 25, 2000
//
// Added the assignment operator.                        Sept. 4, 2000
// Added initialize version of the copy constructor      May.  2, 2001
//
// Fixed assignment operator so that indexing 
// information is not overwritten                        April 22, 2003
//
// Updated stream include to STL                          Jan.  23,2004
// Added initialize()                                     May   30,2005
// Fixed constructor (T*d, long m)                        June  12,2007
//#####################################################################
//
#ifndef __ArrayStructure1D__
#define __ArrayStructure1D__

#ifdef  _DEBUG 
#include <iostream>
#include <cstdio>
using namespace std;
#endif
template <class T>
class ArrayStructure1D 
{
 
    public :

    ArrayStructure1D()
    {
    dataPtr       = 0;
    internalAlloc = 0;
    index1Size    = 0;
    index1Begin   = 0;
    index1End     = 0;
    }

    ArrayStructure1D(T* d, long m)
    {
    dataPtr       = 0;
    internalAlloc = 0;
    initialize(d,m);
    };

    ArrayStructure1D(long m)
    {
    dataPtr       = 0;// will be reset in initialize
    internalAlloc = 0; 
    initialize(m);
    };

    ArrayStructure1D(const ArrayStructure1D& D)
    {
    if(D.dataPtr == 0)
    {
    dataPtr       = 0;
    internalAlloc = 0;
    index1Size    = 0;
    index1Begin   = 0;
    index1End     = 0;
    return;
    }

    index1Size    = D.index1Size;
    index1Begin   = D.index1Begin;
    index1End     = D.index1End;
    dataPtr       = new T[index1Size];
    internalAlloc = 1; 
	long i;
	for(i = 0; i < index1Size; i++){dataPtr[i] = D.dataPtr[i];}
    };

    virtual ~ArrayStructure1D()
    {
    if(internalAlloc == 1) delete [] dataPtr;
    }

    virtual void initialize()
    {
    if(internalAlloc == 1) delete [] dataPtr;
    dataPtr       = 0;
    internalAlloc = 0;
    index1Size    = 0;
    index1Begin   = 0;
    index1End     = 0;
    }

    virtual void initialize(long m)
    {

    if(internalAlloc == 1) 
    {
      if(index1Size != m)
      {
       delete [] dataPtr;
       dataPtr = new T[m];
      }
    }
    else 
    {
      if(dataPtr == 0)
      {
      dataPtr        = new T[m];
      internalAlloc  = 1;
      }
    }
    index1Size    = m;
    index1Begin   = 0;
    index1End     = index1Begin + (index1Size - 1);
    };

    virtual void initialize(T* d, long m)
    {
    if(internalAlloc == 1) delete [] dataPtr;
    internalAlloc = 0;
    dataPtr       = d;
    index1Size    = m;
    index1Begin   = 0;
    index1End     = index1Begin + (index1Size - 1);
    };

    virtual void initialize(const ArrayStructure1D& D)
    {
    if(internalAlloc == 1) 
    {
      if(index1Size != D.index1Size)
      {
       delete [] dataPtr;
       dataPtr = new T[D.index1Size];
      }
    }
    else 
    {
      if(dataPtr == 0)
      {
      dataPtr        = new T[D.index1Size];
      internalAlloc  = 1;
      }
    }
    index1Size    = D.index1Size;
    index1Begin   = D.index1Begin;
    index1End     = D.index1End;
	long i;
	for(i = 0; i < index1Size; i++){dataPtr[i] = D.dataPtr[i];}
    };

    public :

#ifdef _DEBUG 
    T&  operator()(long i1)
    {
    boundsCheck(i1, index1Begin, index1End,1);
    return *(dataPtr +  (i1 - index1Begin));
    };

    const T&  operator()(long i1) const
    {
    boundsCheck(i1, index1Begin, index1End,1);
    return *(dataPtr +  (i1 - index1Begin));
    };

#else

    inline T&  operator()(long i1)
    {
    return *(dataPtr + (i1 - index1Begin));
    };

    inline const T&  operator()(long i1) const
    {
    return *(dataPtr +   (i1 - index1Begin));
    };

#endif

    public :

    T* getDataPointer(){return dataPtr;};

    const T* getDataPointer() const {return dataPtr;};

    void setIndex1Begin(long i) 
    {index1Begin = i; index1End   = index1Begin + (index1Size - 1);};

    long getIndex1Begin() const {return index1Begin;}
    long getIndex1End()   const {return index1End;}
    long getIndex1Size()  const {return index1Size;}

//
//  Get/Set specifically for one dimensional arrays
//
    void setIndexBegin(long i) 
    {index1Begin = i; index1End   = index1Begin + (index1Size - 1);};

    long getIndexBegin() const {return index1Begin;}
    long getIndexEnd()   const {return index1End;}
    long getSize()       const {return index1Size;}

	//
    // Resizes array to exactly newSize
    //
    void resize(long newSize)
	{
	long i;
    T*  newDataPtr = new T[newSize];
	T*  tmpDataPtr;

	if(newSize > index1Size) 
	{
		for(i = 0; i < index1Size; i++) newDataPtr[i] = dataPtr[i];
	}
	else
	{
		for(i = 0; i < newSize; i++) newDataPtr[i] = dataPtr[i];
	}

	index1Size = newSize;
	tmpDataPtr = dataPtr;
	dataPtr    = newDataPtr;

	if(internalAlloc == 1) delete [] tmpDataPtr;
	internalAlloc = 1;

	index1End   = index1Begin + (index1Size - 1);
	}

    void operator=(const ArrayStructure1D& D)
    {
	#ifdef _DEBUG 
	if(index1Size != 0)
	{
    sizeCheck(this->index1Size,D.index1Size);
	}
	#endif
//
//  If null instance, then initialize and acquire indexing
//  from right hand side.
//
    if(index1Size == 0)
    {
    initialize(D.index1Size);
    index1Size    = D.index1Size;
    index1Begin   = D.index1Begin;
    index1End     = D.index1End;
    }

    long i;
    for(i = 0; i < D.index1Size; i++)
    {dataPtr[i] = D.dataPtr[i];}
    }

    void setToValue(T d)
	{
    long i;
	for(i = 0; i < index1Size; i++)
	{dataPtr[i] =  d;}
	}


    void addValue(T d)
	{
    long i;
	for(i = 0; i < index1Size; i++)
	{dataPtr[i] += d;}
	}


    public :

    T*      dataPtr;          // data pointer
    long    index1Begin;      // coordinate 1 starting index
    long    index1End;        // coordinate 1 ending index
    long    index1Size;       // coordinate 1 size
    int     internalAlloc;    // internal allocation flag 

#ifdef _DEBUG 
    static void boundsCheck(long i, long begin, long end, int coordinate)
    {
    if((i < begin)||(i  > end))
    {
    printf("Array index %d out of bounds \n",coordinate);
    printf("Offending index value %ld : Acceptable Range [%ld, %ld] \n",i, begin, end);
    }
    }
#else
static void boundsCheck(long, long, long, int){}
#endif

#ifdef _DEBUG 
    static void sizeCheck(long size1, long size2)
    {
    if(size1 != size2)
    {
    printf("Array Sizes Are Incompatable %ld != %ld \n", size1, size2);
    }
    }
#else
static void sizeCheck(long, long){}
#endif
 
};

#define DoubleArrayStructure1D ArrayStructure1D<double>
#define IntArrayStructure1D    ArrayStructure1D<int>
#define FloatArrayStructure1D  ArrayStructure1D<float>
#define LongArrayStructure1D   ArrayStructure1D<long>
#define CharArrayStructure1D   ArrayStructure1D<char>

#endif 
