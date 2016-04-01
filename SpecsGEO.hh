#ifndef SPECSGEO_hh
#define SPECSGEO_hh
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
class SpecsGEO {
private:
public:
    SpecsGEO();
    virtual ~SpecsGEO();
    virtual double   GetPhi(const double x, const double y,const double philo);
    virtual double   GetPhi(const double phi0,const double philo);
    inline double    GetPhi(const TLorentzVector& tlv,const double philo) { return GetPhi(tlv.Phi(),philo); }
    inline double    GetPhiCLAS(const double x,const double y) { return GetPhi(x,y,-30*TMath::DegToRad()); };
    inline double    GetPhiCLAS(const double phi0) { return GetPhi(phi0,-30*TMath::DegToRad()); };
    inline double    GetPhiRTPC(const double phi0) { return GetPhi(phi0,-90*TMath::DegToRad()); };
    inline double    GetPhiCLAS(const TLorentzVector& tlv) { return GetPhiCLAS(tlv.Phi()); };
    inline double    GetPhiRTPC(const TLorentzVector& tlv) { return GetPhiRTPC(tlv.Phi()); };
    virtual int      GetSector(double fi);
    virtual int      GetSector(const double xx,const double yy);
    inline  int      GetSector(const TLorentzVector& v){ return GetSector(v.Px(),v.Py()); }
    virtual double   DihedralAngle(const TVector3& v1,const TVector3& v2,const TVector3& v3);
    virtual double   DihedralAngle(const TLorentzVector& v1,const TLorentzVector& v2,const TLorentzVector& v3);
    virtual double   CoplanarityAngle(const TVector3& v1,const TVector3& v2,const TVector3& v3);
    virtual double   CoplanarityAngle(const TLorentzVector& v1,const TLorentzVector& v2,const TLorentzVector& v3);
    virtual double   DeltaPhi(const TVector3& v1,const TVector3& v2,const TVector3& v3);
};
SpecsGEO::SpecsGEO(){}
SpecsGEO::~SpecsGEO(){}
double SpecsGEO::GetPhi(const double phi0,const double philo)
{   // return phi in radians in the range [philo,philo+2pi)
    static const double pi=TMath::Pi();
    double phi2=phi0;
    if (phi2<philo) while (phi2<philo)        phi2+=2*pi;
    else            while (phi2>=philo+2*pi)  phi2-=2*pi;
    return phi2;
}
double SpecsGEO::GetPhi(const double xx, const double yy,const double philo)
{   // return phi in radians in the range [philo,philo+2pi)
    // WHY NOT JUST USE ATAN2???
    static const double pi=TMath::Pi();
    double phi0=0;
    if (xx==0) {
	if      (yy==0) phi0 =  0;
        else if (yy>0)  phi0 =  pi/2;
        else            phi0 = -pi/2;
    } else {
        phi0=atan(yy/xx);
        if      (yy>=0 && xx<0) phi0 += pi;
        else if (yy<0  && xx<0) phi0 -= pi;
    }
    return GetPhi(phi0,philo);
}
int SpecsGEO::GetSector(double fi)
{
    // returns sector in range 0-5 (user must add 1 to get CLAS sectors)
    static const double pi=TMath::Pi();
    while (fi<0) fi += 2*pi;
    return (int)((fi+pi/6)/(pi/3))%6;
}
int SpecsGEO::GetSector(const double xx,const double yy)
{
    return GetSector(GetPhi(xx,yy,-TMath::Pi()/6));
}
double SpecsGEO::DihedralAngle(const TVector3& v1,const TVector3& v2,const TVector3& v3)
{
    return atan2 (  v2.Mag()*v1.Dot(v2.Cross(v3)) ,
                 (v1.Cross(v2)).Dot(v2.Cross(v3)) );
}
double SpecsGEO::DihedralAngle(const TLorentzVector& v1,const TLorentzVector& v2,const TLorentzVector& v3)
{
    return DihedralAngle(v1.Vect(),v2.Vect(),v3.Vect());
}
double SpecsGEO::CoplanarityAngle(const TVector3& v1,const TVector3& v2,const TVector3& v3)
{
    const TVector3 n12=v1.Cross(v2).Unit();
    const TVector3 n23=v2.Cross(v3).Unit();
    return acos(n12.Dot(n23));
}
double SpecsGEO::CoplanarityAngle(const TLorentzVector& v1,const TLorentzVector& v2,const TLorentzVector& v3)
{
    return CoplanarityAngle(v1.Vect(),v2.Vect(),v3.Vect());
}
double SpecsGEO::DeltaPhi(const TVector3& v1,const TVector3& v2,const TVector3& v3)
{
    const double yy = v1.Dot( v2.Cross( v1.Cross(v3) ) );
    const double xx = v1.Mag() * v2.Dot( v1.Cross(v3) );
    return GetPhi(xx,yy,0);
}
#endif
