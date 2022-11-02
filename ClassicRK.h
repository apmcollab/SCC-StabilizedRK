
//
//####################################################################
//                    ClassicRK.h 
//####################################################################
/**
  
   A templated class for Classic Runge-Kutta methods of 
   order 1 - 4.<p>
  
   The minimal functionality required of the classes 
   that are used in this template
   <pre> 
   State :
  
   State()                           (null constructor)
   State(const State&)               (copy constructor)
  
   operator =                        (duplicate assignemnt)
   operator +=                       (incremental addition)
   operator *(double alpha)          (scalar multiplication)
  
   StateOperator :

   void apply(State&S, State &FS)
   </pre>
<i>Source</i>: 
<A HREF="../ClassicRK.h">ClassicRK.h</A><p>

@author Chris Anderson (C) UCLA 
@version May 18, 2000
*/
//#####################################################################
// Chris Anderson (C) UCLA                               Feb. 23, 2006
//#####################################################################
//

#ifndef CLASSIC_RK_
#define CLASSIC_RK_

template <class State, class StateOperator> class ClassicRK
{
public :
    
ClassicRK(StateOperator& F, State& S, long Order) 
: Sn(S),Snp1(S),Sstar(S),FS(S)
{
    Fop     = &F;
    rkOrder = Order;
    setRKcoefficients();
}

void initialize(StateOperator& F, State& S, long Order) 
{
    Fop     = &F;
    rkOrder = Order;
    setRKcoefficients();

    // initialize temporaries

    Sn.initialize(S);
    Snp1.initialize(S);
    Sstar.initialize(S);
    FS.initialize(S);
}

ClassicRK() 
: Sn(),Snp1(),Sstar(),FS()
{
    Fop     = 0;
    rkOrder = 1;
    setRKcoefficients();
}

void setOperator(StateOperator& F)
{Fop = &F;}

void setOrder(long order)
{
    rkOrder = order; 
    setRKcoefficients();
};

long getOrder() {return rkOrder;};

void advance(State& YnIn, State& FYnIn,double dt,State& YnOut, State& FYnOut)
{
    long i;

    Sn     = YnIn;
    Snp1   = YnIn;
    Sstar  = YnIn;
    FS     = FYnIn; 
//
// Loop  over stages
//
    for (i=0; i < rkOrder; i++)
    {
        if(i > 1)
        {
        Fop->apply(Sstar,FS);
        }

        Snp1.axpy(dt*rkk[i],FS);     // Snp1 += FS*(dt*rkk[i]);

        if((i+1) < rkOrder)
        {
          Sstar   = Sn;
          Sstar.axpy(dt*rk[i+1],FS); // Sstar  += FS*(dt*rk[i+1]);
        }
    }

    YnOut = Snp1;
    Fop->apply(YnOut,FYnOut);
}


void advance(State& S, double dt)
{
    long i;

    Sn     = S;
    Snp1   = S;
    Sstar  = S;
    FS     = S; // Assignment to insure that FS is initialized
//
// Loop  over stages
//
    for (i=0; i < rkOrder; i++)
    {
        Fop->apply(Sstar,FS);

        Snp1.axpy(dt*rkk[i],FS);     // Snp1 += FS*(dt*rkk[i]);

        if((i+1) < rkOrder)
        {
          Sstar   = Sn;
          Sstar.axpy(dt*rk[i+1],FS); // Sstar  += FS*(dt*rk[i+1]);
        }
    }

    S = Snp1;
}

private :

    State Sn;
    State Snp1;
    State Sstar;
    State FS;

    StateOperator* Fop;

    long  rkOrder;
    double  rk[4];
    double rkk[4];
//
//  Utility Functions
//
void setRKcoefficients()
{
    if(rkOrder < 1) {rkOrder = 1;}
    if(rkOrder > 4) {rkOrder = 4;}

    switch(rkOrder)
    {
        case 1:  rk[0]  = 0.0;

                 rkk[0] = 1.0;
        break;
        case 2:  rk[0]  = 0.0;     rk[1]  = 1.0;

                 rkk[0] = 0.5;     rkk[1] = 0.5;
        break;
        case 3:  rk[0]  = 0.0;     rk[1]  =  1.0/3.0;  
                 rk[2]  = 2.0/3.0;

                 rkk[0] =.25;      rkk[1] =  0.0;      
                 rkk[2] =.75;
        break;
        case 4:  rk[0]  = 0.0;      rk[1]  = 1.0/2.0;   
                 rk[2]  = 1.0/2.0;  rk[3] = 1.0;

                rkk[0]  = 1.0/6.0; rkk[1] = 1.0/3.0;  
                rkk[2]  = 1.0/3.0; rkk[3] = 1.0/6.0;
        break;
    } 
}
};
#endif
 
