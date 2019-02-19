/* 
Evaluator classes and tests

Goals:

- given a function with options
    f(x, opt1, opt2, ...)
  encapsule it into an object that works as
    F f(opt1, opt2, ...);
    f(x);

- given functions
    f1(x), f2(x)
  define operators that gives the functional sum, product... 
  which works as
    F f = f1 * f2;
    f(x);

This file only serves as an example
Use g++ for compilation
*/


#include <iostream>


/* Evaluator Base Class */

class Evaluator {
public:
    virtual float evaluate(float x)=0;
    float operator() (float x) { 
        return evaluate(x); 
    }
};


/* Evaluator Mixers */

class EvaluatorBinaryMixer: public Evaluator {
public:
    EvaluatorBinaryMixer(): 
        _p1(nullptr), _p2(nullptr) {}
    EvaluatorBinaryMixer(Evaluator& e1, Evaluator& e2): 
        _p1(&e1), _p2(&e2) {}
    EvaluatorBinaryMixer(const EvaluatorBinaryMixer& e):
        _p1(e._p1), _p2(e._p2) {}
    EvaluatorBinaryMixer& operator=(const EvaluatorBinaryMixer& e) {
        _p1 = e._p1;
        _p2 = e._p2;
        return *this;
    }
    virtual ~EvaluatorBinaryMixer() {}
protected:
    Evaluator* _p1;
    Evaluator* _p2;
};

class EvaluatorProduct: public EvaluatorBinaryMixer {
public:
    EvaluatorProduct(): 
        EvaluatorBinaryMixer() {}
    EvaluatorProduct(Evaluator& e1, Evaluator& e2): 
        EvaluatorBinaryMixer(e1, e2) {}
    EvaluatorProduct(const EvaluatorProduct& e):
        EvaluatorBinaryMixer(e) {}
    EvaluatorProduct& operator=(const EvaluatorProduct& e) {
        EvaluatorBinaryMixer::operator=(e);
        return *this;
    }
    virtual float evaluate(float x) {
        return _p1->evaluate(x) * _p2->evaluate(x);
    }
};

EvaluatorProduct operator* (Evaluator& e1, Evaluator& e2) {
    EvaluatorProduct p(e1, e2);
    return p;
}



/* Test Case */

class LinearNodal: public Evaluator {
public:
    LinearNodal(float x_mid, float x_min, float x_max): 
        _mid(x_mid), _min(x_min), _max(x_max) {}
    LinearNodal(int A, const float* mesh): 
        _mid(mesh[A]), _min(mesh[A-1]), _max(mesh[A+1]) {}
    virtual float evaluate(float x) {
        if (x < _min || x > _max)
            return 0.;
        if (x < _mid)
            return (x - _min) / (_mid - _min);
        else 
            return (_max - x) / (_max - _mid);
    }
private:
    float _mid, _min, _max;
};


int main() {

    float mesh[4] = {0.,1.,2.,3.};

    LinearNodal N1(1, mesh);
    LinearNodal N2(2, mesh);
    EvaluatorProduct p = N1 * N2;
    // Evaluator* p = new EvaluatorProduct(N1, N2);

    for (int i=0; i<10; i++) {
        float x = i * 3./10.;
        std::cout << x << ": "
            << N1(x) << ", "
            << N2(x) << ", "
            << p(x) << std::endl;
            // << p->evaluate(x) << std::endl;
    }
    // delete p;

    return 0;
}

