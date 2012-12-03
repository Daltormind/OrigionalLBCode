#include "wet.h"

void wet::Dimensionlessnum(void)
{

    cout << " dt,dx,mo " << dt << " " << dx << " " << mo << endl;
    
    vl=(tauliquid-.5)/3;
    vg=(taugas-.5)/3;
    
    cout << "phi11=" << phi11 << endl;
    
    w=sqrt(2.0*kappa_p/A)*(-phi11);
    
    cout << "w=" << w << endl;
    
    stl=sqrt(2*kappa_p*A/9)*(1-pow(1+w,1.5));
    
    cout << "stl=" << stl << endl;
    
    stg=sqrt(2*kappa_p*A/9)*(1-pow(1-w,1.5));
    
    cout << "stg=" << stg <<endl;
    
    we=1.0*(pow(initUX,2)+pow(initUY,2)+pow(initUZ,2))*2*dropletR*dx/sqrt(8*kappa_p*A/9);
    
    cout << "we=" << we << endl;
    
    Re=sqrt(pow(initUX,2)+pow(initUY,2)+pow(initUZ,2))*2*dropletR/vl;
    
    cout << "Re=" << Re << endl;
    Oh=vl*sqrt(1.0)/sqrt(sqrt(8*kappa_p*A/9)*2*dropletR);
    
    cout << "Oh=" << Oh << endl;
    
    Mo=gama*(tau1-0.5);
    
    Pe=sqrt(kappa_p)*sqrt(pow(initUX,2)+pow(initUY,2)+pow(initUZ,2))/Mo/pow(A,1.5);
    
    cout << "Pe=" << Pe << endl;
    







}
