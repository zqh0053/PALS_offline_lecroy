#ifndef GETPARAMETERS_H
#define GETPARAMETERS_H 1
#include "runaction.h"
class getparameters
{
public:
    getparameters(RunAction *);
    ~getparameters();
    void readfile(string );
private:
    RunAction *frunaction;
};


#endif
