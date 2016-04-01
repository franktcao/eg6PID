#ifndef SPECSFID_hh
#define SPECSFID_hh
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TCutG.h"
#include "SpecsGEO.hh"
class SpecsFID {
private:
public:
    SpecsGEO GEO;
    SpecsFID();
    virtual ~SpecsFID();
    virtual bool BostedIC(const double x,const double y,const int ip,const bool loose);
    virtual bool BostedIC(const int irun,const double q,const double p,const double cx,const double cy,const double vz,const bool loose);
    virtual bool BostedIC(const int irun,const double q,const TLorentzVector& v,const double vz,const bool loose);
    virtual bool SOL(const double vz,const double cz);
    virtual bool RTPCnose(double vz,double cz);
    virtual bool RTPCbarrel(double vz,double cz);
    virtual bool RTPCendplate(double vz,double cz);
    virtual bool RTPCring(double vz,double cz);
    virtual TVector3 ECxyz2uvw(const TVector3& xyz);
    virtual bool EC(const double x,const double y,const double z);
    virtual bool MohammadIC(const double xx,const double yy);
    virtual bool MohammadIC(const TVector3& r1dir,const TVector3& r1pos);
    virtual TVector3 ProjectDC1toIC(const TVector3& r1dir,const TVector3& r1pos,const double z0);
    virtual bool FxIC(const double xx,const double yy);
    virtual bool FxInSector(const double X,const double Y);
    virtual bool FxInSector(const TVector3& r1dir,const TVector3& r1pos);
    virtual bool RejectHotIC(const double xx,const double yy);
    virtual TVector3 Bosted2IC(const TLorentzVector& v,const double vz,const double q);
};
SpecsFID::SpecsFID(){}
SpecsFID::~SpecsFID(){}

bool SpecsFID::BostedIC(const double x,const double y,const int ip,const bool loose)
{
    // Translated from:  http://www.jlab.org//Hall-B/secure/eg1-dvcs/software/fid_cut_ic.f
    // Modification:
    // Input now requires z-vertex relative to IC. EG6 needs that due to long target.
    // Need to check solenoid polarity and ipart for EG6.

    // return value: wether track intersected IC

    // inputs:
    // x,y is hit position on IC (cm)
    // ip is 0,1,2 for part A,B,C
    // loose is true for regular cut, false for tight cut

    // these  define  4 of the  octogon  faces  (vertical={horizontal)
    static const double xp[3] ={ 17.0,  18.5,  18.5};
    static const double xm[3] ={-17.0, -17.5, -17.5};
    static const double yp[3] ={ 18.0,  18.0,  18.0};
    static const double ym[3] ={-20.0, -20.0, -20.0};
    // these  define  the  diagonal  sides
    static const double xdp[3] ={ 26.0,  26.0,  26.0};
    static const double xdm[3] ={-27.0, -27.0, -27.0};
    static const double ydp[3] ={ 23.0,  23.0,  23.0};
    static const double ydm[3] ={-25.0, -23.0, -23.0};
    // these  define  the  bottom  plate
    static const double xbm[3] ={0., -8.0, -8.5};
    static const double xbp[3] ={0.,  8.0,  8.5};

    // tight or loose cut, extra cm for tigth
    const double off = loose ? 0.0 : 1.0;

    // inside octagon
    const bool ic1 = (x<xp[ip]+off) && (x>xm[ip]-off) &&
        (y<yp[ip]+off) && (y>ym[ip]-off) &&
        x < (xdp[ip]+off)*(1.-y/(ydp[ip]+off)) &&
        x < (xdp[ip]+off)*(1.-y/(ydm[ip]-off)) &&
        x > (xdm[ip]-off)*(1.-y/(ydp[ip]+off)) &&
        x > (xdm[ip]-off)*(1.-y/(ydm[ip]-off)) ;

    // inside thick bottom plate
    const bool ic2 = y>ym[ip]-off  && y<0 && 
        x>xbm[ip]-off && x<xbp[ip]+off;

    return ic1 || ic2;
}

bool SpecsFID::BostedIC(const int irun,
                          const double q,
                          const double p,
                          const double cx,const double cy,
                          const double vz,
                          const bool loose)
{
    // Translated from:  http://www.jlab.org//Hall-B/secure/eg1-dvcs/software/fid_cut_ic.f
    // Modification:
    // Input now requires z-vertex relative to IC. EG6 needs that due to long target.
    // Need to check solenoid polarity and ipart for EG6.
    
    // return value: wether track intersected IC

    // inputs:
    // irun = run number
    // q = charge (-1,0,1)
    // p = momentum in GeV/c
    // cx,cy = direction cosines at target vertex
    // vz = vz-zIC (cm)  -- for EG6, this is just vz since IC is at 0
    // regular = (0=tight, 1=loose)

    // need to check this for eg6 (solenoid polarity?):
    const double targsign= -1;
    //const double targsign= irun>60222 ? 1 : -1;
    
    // need to check this for eg6:
    int ipart=0;
    if      (irun>=59221) ipart=1;
    else if (irun>=60220) ipart=2;

    //const double dphi = q*targsign*0.1845/p; // where 0.1845 from?
    const double dphi = q*targsign*0.1545/p;

    const double x0 = fabs(vz)*cx;
    const double y0 = fabs(vz)*cy;

    const double xic = x0*cos(dphi)+y0*sin(dphi);
    const double yic = y0*cos(dphi)-x0*sin(dphi);

    return BostedIC(xic,yic,ipart,loose);
}
TVector3 SpecsFID::Bosted2IC(const TLorentzVector& v,const double vz,const double q)
{
    const double p=v.P();
    const double cx=v.Px()/p;
    const double cy=v.Py()/p;

    const double targsign=-1;
    const double dphi = q*targsign*0.1845/p; // where 0.1845 from?

    const double x0 = fabs(vz)*cx;
    const double y0 = fabs(vz)*cy;

    const double xic = x0*cos(dphi)+y0*sin(dphi);
    const double yic = y0*cos(dphi)-x0*sin(dphi);

    return TVector3(xic,yic,0);
}
bool SpecsFID::BostedIC(const int irun,
                          const double q,
                          const TLorentzVector& v,
                          const double vz,
                          const bool loose)
{
    return BostedIC(irun,q,v.P(),v.Px()/v.P(),v.Py()/v.P(),vz,loose);
}
bool SpecsFID::MohammadIC(const double xx,const double yy)
{
  static const int nn=11;
  static const double xpos[nn]={ -11.15, -11.15, -23.1, -23.1, -10.3,   9.91, 23.73, 23.73, 12.3,   12.3,  -11.5};
  static const double ypos[nn]={ -26.07, -23.1,  -12.85, 11.5,  22.95, 22.95, 13.1, -12.4, -22.36, -26.07, -26.07};
  static TCutG* geocut=new TCutG("geocut", nn-1, xpos,ypos);
  return !geocut->IsInside(xx,yy);
}
TVector3 SpecsFID::ProjectDC1toIC(const TVector3& r1dir,const TVector3& r1pos,const double z0)
{
    // this projects to a x/y plane normal to z
    const TVector3 shift = (r1pos.Z()-z0) / r1dir.Z() * r1dir;
    return r1pos-shift;
}
bool SpecsFID::MohammadIC(const TVector3& r1dir,const TVector3& r1pos)
{
    const TVector3 icpos=ProjectDC1toIC(r1dir,r1pos,16.0);  // <-- Mohammad's cut is at the downstream end of IC(?)
    return MohammadIC(icpos.X(),icpos.Y());
}
bool SpecsFID::FxInSector(const double X,const double Y)
{
    const int sec=GEO.GetSector(X,Y)+1;
    if (sec%6<3)
    {
        if (Y<X*tan(TMath::Pi()*((sec-1)/3.-1./9))) return 0;
        if (Y>X*tan(TMath::Pi()*((sec-1)/3.+1./9))) return 0;
    }
    else
    {
        if (Y>X*tan(TMath::Pi()*((sec-1)/3.-1./9))) return 0;
        if (Y<X*tan(TMath::Pi()*((sec-1)/3.+1./9))) return 0;
    }
    return 1;
}
bool SpecsFID::FxInSector(const TVector3& r1dir,const TVector3& r1pos)
{
    const TVector3 icpos=ProjectDC1toIC(r1dir,r1pos,16.0);
    return FxInSector(icpos.X(),icpos.Y());
}
bool SpecsFID::FxIC(const double xx,const double yy)
{
    // inputs are xc,yc from ICPB
    // this is to reject gammas near the inner/outer edges of the IC
    static const double dx=1.346; // cm
    static const double dy=1.360; // cm
    static const double nin=3.25;
    static const double nout=10.75;
    static const double root2=sqrt(2);

    // INNER:
    if (fabs(xx)/dx <= nin  &&
        fabs(yy)/dy <= nin  &&
        fabs(xx/dx - yy/dy) <= nin*root2 &&
        fabs(xx/dx + yy/dy) <= nin*root2 ) return 0;

    // OUTER:
    if (fabs(xx)/dx >= nout  ||
        fabs(yy)/dy >= nout  ||
        fabs(xx/dx - yy/dy) >= nout*root2 ||
        fabs(xx/dx + yy/dy) >= nout*root2 ) return 0;
    
    return 1;
}
bool SpecsFID::SOL(const double vz,const double cz)
{
    // returns wether track cleared the solenoid structure
    static const double zsol = -64.0 + 20.96 / 2.0;
    return cz > cos(atan2(11., zsol-vz));
}


bool SpecsFID::RTPCnose(double vz,double cz)
{
    // returns wether track hits the upstream target holder nose
    static const double radius=2.5; // radius of target (mm)
    static const double zpos=-84; // downstream end of support tube (mm)
    vz = vz*10+640; // convert to RTPC coordinates
    return cz < cos( atan2( radius, zpos-vz ) );
}
bool SpecsFID::RTPCbarrel(double vz,double cz)
{
    // returns wether track hits the barrel region
    static const double radius=72.65;
    vz = vz*10+640; // convert to RTPC coordinates
    return cz < cos( atan2( radius, 100-vz ) );
}
bool SpecsFID::RTPCendplate(double vz,double cz)
{
    // returns wether track hits downstream endplate
    static const double radius=60.09;
    vz = vz*10+640; // convert to RTPC coordinates
    return cz > cos( atan2( radius, 100-vz ) );
}
bool SpecsFID::RTPCring(double vz,double cz)
{
    return !RTPCbarrel(vz,cz) && !RTPCendplate(vz,cz);
}


TVector3 SpecsFID::ECxyz2uvw(const TVector3& xyz)
{
    static const double zoffset =  510.32;
    static const double ec_the  =    0.43633230;
    static const double ylow    = -182.97400000;
    static const double yhi     =  189.95600000;
    static const double tgrho   =    1.95325000;//1.097620829
    static const double sinrho  =    0.89012560;
    static const double cosrho  =    0.45571500;
    static const double sinthe  = sin(ec_the);
    static const double costhe  = cos(ec_the);

    double clas_phi = xyz.Phi()*180/TMath::Pi();
    if(clas_phi < 0.) clas_phi += 360;
    clas_phi += 30;
    if(clas_phi >= 360.) clas_phi -= 360;
    const double sector_phi = (int)(clas_phi/60.)*(TMath::Pi()/3);

    double rot[3][3];
    rot[0][0] =  costhe*cos(sector_phi);
    rot[0][1] = -sin(sector_phi);
    rot[0][2] =  sinthe*cos(sector_phi);
    rot[1][0] =  costhe*sin(sector_phi);
    rot[1][1] =  cos(sector_phi);
    rot[1][2] =  sinthe*sin(sector_phi);
    rot[2][0] = -sinthe;
    rot[2][1] =  0.;
    rot[2][2] =  costhe;

    double xyzi[3];
    xyzi[1] = xyz.X()*rot[0][0]+xyz.Y()*rot[1][0]+xyz.Z()*rot[2][0];
    xyzi[0] = xyz.X()*rot[0][1]+xyz.Y()*rot[1][1]+xyz.Z()*rot[2][1];
    xyzi[2] = xyz.X()*rot[0][2]+xyz.Y()*rot[1][2]+xyz.Z()*rot[2][2] - zoffset;

    TVector3 uvw;
    uvw.SetX( (xyzi[1] - ylow)/sinrho );
    uvw.SetY( (yhi-ylow)/tgrho - xyzi[0] + (yhi-xyzi[1])/tgrho );
    uvw.SetZ( ( (yhi-ylow)/tgrho + xyzi[0] + (yhi-xyzi[1])/tgrho ) / 2. / cosrho );

    return uvw;
}
bool SpecsFID::EC(const double x,const double y,const double z)
{
    const TVector3 uvw=ECxyz2uvw(TVector3(x,y,z));
    if (uvw.X()<40 || uvw.X()>400) return 0;
    if (uvw.Y()>360) return 0;
    if (uvw.Z()>390) return 0;
    return 1; 
}
bool SpecsFID::RejectHotIC(const double xgl,const double ygl)
{
    static const double dx=1.346; // cm
    static const double dy=1.360; // cm
    static const int nhot=9;
    static const int pix[nhot*2]=
    {
         0, 0,
         3, 4,
        -8,-2,
        -4,-6,
        -2,-6,
        -1,-8,
         3,-10,
        -5, 8,
        -9,-6
    };
    
    const int xpix=(int)round(xgl/dx);
    const int ypix=(int)round(ygl/dy);

    for (int ii=0; ii<nhot; ii++)
    {
        if (xpix==pix[ii*2] && ypix==pix[ii*2+1]) return 0;
    }
    return 1;
}
#endif
