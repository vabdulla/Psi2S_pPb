#ifndef bin_h
#define bin_h

#include <utility>
#include <tuple>
#include <iostream>
#include <set>

using namespace std;

// a simple template class to store a bin and overload the equality operator

// define a few common uses of the template class
template <typename T> class bin : public pair<T,T> {
   public:
      bin(T a, T b) : pair<T,T>(b,a) {};
      T low() const {return this->second;}
      T high() const {return this->first;}
};
typedef bin<double> binD;
typedef bin<float> binF;
typedef bin<int> binI;

// associate three such bins to make an analysis bin
class anabin : public tuple<binF,binF,binI> {
   public:
      anabin(float rapmin, float rapmax, float ptmin, float ptmax, int centmin, int centmax) :
         tuple<binF,binF,binI> (binF(rapmin,rapmax), binF(ptmin, ptmax), binI(centmin, centmax)) {};
      binF rapbin() const {return get<0>(*this);};
      binF ptbin() const {return get<1>(*this);};
      binI centbin() const {return get<2>(*this);};
      void setrapbin(binF rapbin) {get<0>(*this) = rapbin;};
      void setptbin(binF ptbin) {get<1>(*this) = ptbin;};
      void setcentbin(binI centbin) {get<2>(*this) = centbin;};
      void print() const {
         cout << "rap=[" << get<0>(*this).low() << "," << get<0>(*this).high() <<
            "], pt=[" << get<1>(*this).low() << "," << get<1>(*this).high() <<
            "], cent=[" << get<2>(*this).low() << "," << get<2>(*this).high() << "]" << endl;
      }
};

set<anabin> allbins() {
   set<anabin> ans;
   // HIN-16-004
   //cout<<"rapmin: "<<rap->getMin()<<"rapmax: "<<rap->getMax()<<"ptmin  "<<pt->getMin()<<"ptmax  "<<pt->getMax()<<endl;
   ans.insert(anabin(1.46,1.93,3.0,4.0,0,200));
   ans.insert(anabin(1.46,1.93,4.0,5.0,0,200));
   ans.insert(anabin(1.46,1.93,5.0,6.5,0,200));
   ans.insert(anabin(1.46,1.93,4.0,6.5,0,200));
   ans.insert(anabin(1.46,1.93,6.5,10.0,0,200));
   ans.insert(anabin(1.46,1.93,10.0,30.0,0,200));
   ans.insert(anabin(1.46,1.93,3.0,6.5,0,200));
   ans.insert(anabin(1.46,1.93,4.0,10.0,0,200));
   ans.insert(anabin(1.46,1.93,4.0,30.0,0,200));
   ans.insert(anabin(1.46,1.93,6.5,30.0,0,200));

   ans.insert(anabin(1.03,1.46,5.0,6.5,0,200));
   ans.insert(anabin(1.03,1.46,6.5,10.0,0,200));
   ans.insert(anabin(1.03,1.46,5.0,10.0,0,200));
   ans.insert(anabin(1.03,1.46,10.0,30.0,0,200));
   ans.insert(anabin(1.03,1.46,5.0,30.0,0,200));
   ans.insert(anabin(1.03,1.46,6.5,30.0,0,200));

   ans.insert(anabin(0.43,1.03,6.5,10.0,0,200));
   ans.insert(anabin(0.43,1.03,10.0,30.0,0,200));
   ans.insert(anabin(0.43,1.03,6.5,30.0,0,200));

   ans.insert(anabin(-0.47,0.43,6.5,10.0,0,200));
   ans.insert(anabin(-0.47,0.43,10.0,30.0,0,200));
   ans.insert(anabin(-0.47,0.43,6.5,30.0,0,200));

   ans.insert(anabin(-1.37,-0.47,6.5,10.0,0,200));
   ans.insert(anabin(-1.37,-0.47,10.0,30.0,0,200));
   ans.insert(anabin(-1.37,-0.47,6.5,30.0,0,200));
   
   ans.insert(anabin(-1.97,-1.37,6.5,10.0,0,200));
   ans.insert(anabin(-1.97,-1.37,10.0,30.0,0,200));
   ans.insert(anabin(-1.97,-1.37,6.5,30.0,0,200));
   
   ans.insert(anabin(-2.4,-1.97,4.0,5.0,0,200));
   ans.insert(anabin(-2.4,-1.97,5.0,6.5,0,200));
   ans.insert(anabin(-2.4,-1.97,6.5,10.0,0,200));
   ans.insert(anabin(-2.4,-1.97,5.0,10.0,0,200));
   ans.insert(anabin(-2.4,-1.97,10.0,30.0,0,200));
   ans.insert(anabin(-2.4,-1.97,4.0,6.5,0,200));
   ans.insert(anabin(-2.4,-1.97,4.0,10.0,0,200));
   ans.insert(anabin(-2.4,-1.97,6.5,30.0,0,200));
   ans.insert(anabin(-2.4,-1.97,4.0,30.0,0,200));
   
   ans.insert(anabin(-2.07,1.13,6.5,30.0,0,200));
   ans.insert(anabin(-2.07,1.13,6.5,10.0,0,200));
   ans.insert(anabin(-2.07,1.13,10.0,30.0,0,200));    
   
//ans.insert(anabin(1.6,2.0,6.5,50,0,200));
//ans.insert(anabin(2.0,2.4,6.5,50,0,200));

   return ans;
};

#endif // #ifndef bin_h
