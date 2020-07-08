//
// DynamicVectorArray.h 
//
// Author: Chris Anderson  
// (C) UCLA 2012 
//
// DynamicVectorArray.h: 
//
//////////////////////////////////////////////////////////////////////

#ifndef __DynamicVectorArray__
#define __DynamicVectorArray__

#include <stdio.h>
template <class Vector> class DynamicVectorArray  
{
public:

	DynamicVectorArray()
	{
		vArray       =  0;
		arraySize    =  0;
	};

	virtual ~DynamicVectorArray()
	{
		destroyData();
	};

    DynamicVectorArray(long n)
	{
	    vArray    = 0;
		arraySize = 0;
		initialize(n);
	}

    DynamicVectorArray(long n, const Vector& V)
	{
	    vArray    = 0;
		arraySize = 0;
		initialize(n,V);
	}

    DynamicVectorArray(const DynamicVectorArray < Vector > & V)
	{
        vArray    = 0;
		arraySize = 0;
    	if(V.arraySize == 0) {return;}
		arraySize = V.arraySize;
		vArray     = new Vector*[arraySize];
 
		long i;

        for(i = 0; i < arraySize; i++)
        {
        vArray[i] = new Vector(V[i]);
        }
	}

    void expand(long nAdditional)
    {
       long i;
       Vector** vArrayTmp = new Vector*[arraySize + nAdditional];

       if(vArray != 0)
       {
            for(i = 0; i < arraySize; i++) vArrayTmp[i] = vArray[i];
            delete [] vArray;
            vArray = vArrayTmp;
            for(i = arraySize; i < arraySize + nAdditional; i++) vArray[i] = new Vector(*vArray[0]);
            arraySize += nAdditional;
       }
       else
       {
            arraySize = nAdditional;
            vArray    = vArrayTmp;
            for(i = 0; i < arraySize; i++) vArray[i] = new Vector();
       }
    }


    void initialize()
    {
     if(vArray != 0) destroyData();
     vArray       =  0;
     arraySize    =  0;
    }

	void initialize(long n)
	{
        long i;

		if(vArray != 0) destroyData();
		vArray     = new Vector*[n];
        arraySize = n;

        for(i = 0; i < arraySize; i++)
        {
        vArray[i] = new Vector();
        }
	}


 	void initialize(long n,const Vector& V)
	{
	    long i; 

		if(vArray != 0) destroyData();
		vArray     = new Vector*[n];
	    arraySize = n;

		for(i = 0; i < arraySize; i++)
        {vArray[i] = new Vector(V);}
	}

	void initialize(const DynamicVectorArray < Vector > & V)
	{
		if(vArray != 0) destroyData();
		arraySize = V.arraySize;
		vArray     = new Vector*[arraySize];
 
		long i;

        for(i = 0; i < arraySize; i++)
        {
        vArray[i] = new Vector(V[i]);
        }
	}

    long getSize() const {return arraySize;}

    Vector** getDataPtr(){return vArray;}

#ifdef _DEBUG 	
	Vector& operator[](long i)
	{
	boundsCheck(i,arraySize);
	return *vArray[i];
	}

    const Vector& operator[](long i) const
	{
	boundsCheck(i,arraySize);
	return *vArray[i];
	}
#else
	inline Vector& operator[](long i)
	{
	return *vArray[i];
	}

    inline const Vector& operator[](long i) const
	{
	return *vArray[i];
	}
#endif

   void destroyData()
   {
        long i;
		if(vArray != 0) 
		{
        for(i = 0; i < arraySize; i++) delete vArray[i];
        delete [] vArray;
		}
   }

	long arraySize;
	Vector** vArray;

#ifdef _DEBUG 
    static void boundsCheck(long i,long aSize)
    {
    if((i < 0)||(i  >= aSize))
    {
    printf("VectorArray index out of bounds \n");
    printf("Offending index value %ld : Acceptable Range [%d, %ld] \n",i, 0, aSize-1);
    }
    }
#else
static void boundsCheck(long, long){}
#endif


};
#endif  
 
