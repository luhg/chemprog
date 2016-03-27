# This "awk" program is used to calculate the pair distribution function g(r)
# and the "bonded" pair coordination function g0(r) as well as their runing
# coordination numbers n(r) and n0(r) from the xyz file of some MD packages,
# such as CPMD and CP2K.
# Usage:
#       gawk -f rdf.awk name_of_xyz_file
# This program needs a file "pair.txt" in the current directory to specific
# the parameters as you require.
# The format and means is list in below [use gKO(r)as an example]:
# atoma O   the center atom
# atomb H   the coordinated atom
# cellx 11.2827   The size of cell:x
# celly 11.2827   The size of cell:y
# cellz 11.2827   The size of cell:z
# deltar 0.02     The radial interval (in Angstrom)
# ntotal 500      The total radial step 
# stepmin  50000   from which dynamics steps you want to consider
# stepmax 100000   to which dynamics steps you want to consider
# bonded 1       if calcualting the "bonded" pair coordination function: 0: no; 1: yes
# Note: if bonded=1, the atoma and atomb shoud be H and O to extract the H-bonds in water.
BEGIN{
#Set a pair of atoms you want to consider.
     atoma="H";
     atomb="O";
#Set the size of cell: x, y, z.
     cellx=12.4138;
     celly=12.4138;
     cellz=12.4138;
#Set the minimum and maxinum dynamics steps you want to calculated.
     stepmin=10001;
     stepmax=20000;
#Set the interval (delta r, in Angstrom) and the total steps of radius (r).
     deltar=0.05;
     ntotal=200;
     bonded=1;
# Read user's initial parameter from the file "pair.txt" in current directory.
     for (i=1;i<=100;i++)
       {
       getline mypar < "pair.txt";
       split(mypar,mypar0);
       if ( mypar0[1] == "atoma" ) atoma=mypar0[2];
       if ( mypar0[1] == "atomb" ) atomb=mypar0[2];
       if ( mypar0[1] == "cellx" ) cellx=mypar0[2];
       if ( mypar0[1] == "celly" ) celly=mypar0[2];
       if ( mypar0[1] == "cellz" ) cellz=mypar0[2];
       if ( mypar0[1] == "stepmin" ) stepmin=mypar0[2];
       if ( mypar0[1] == "stepmax" ) stepmax=mypar0[2];
       if ( mypar0[1] == "deltar" ) deltar=mypar0[2];
       if ( mypar0[1] == "ntotal" ) ntotal=mypar0[2];
       if ( mypar0[1] == "bonded") bonded=mypar0[2];
       }
# Set the initial value of count array.
     for (i=1;i<=ntotal;i++)
       {
       pair00[i]=0;
       pair11[i]=0;
       }
# The parameter for the first read content in cell. 
     nread=0;
     nstep=0;
     dist=deltar*ntotal;
# Change the order of 3x3x3 cells.
     ncell[1]=0;
     ncell[2]=-1;
     ncell[3]=+1;
print "Now, we will calculate the pair correlation function between ",atoma," and ", atomb,".";
print "The cell size is x=",cellx,", y=",celly,"and z=",cellz,".";
print "The dynamic steps are from ",stepmin,"to",stepmax,".";
print "The inteval of distance (delta r, in Angstrom) is,"deltar",", "and the total step is ",ntotal,".";
print "If the above parameters have some thing wrong, please change your file 'pait.txt'to correct them." 
nconfig=-1;
}
{
if (nconfig == -1){headxyz=$0;natom=$1;nconfig=0;}
if ($0 == headxyz) {
nconfig=nconfig+1;
getline;
# read the coordinates of atoms
na=0;
nb=0;
for (i=1;i<=natom;i++)
   {
   getline;
   atomx[i]=$2;
   atomy[i]=$3;
   atomz[i]=$4;
   if ( $1 == atoma ) 
     {
     na=na+1;
     ax[na]=$2;
     ay[na]=$3;
     az[na]=$4;
      }
   if ( $1 == atomb ) 
     {
     nb=nb+1;
     bx[nb]=$2;
     by[nb]=$3;
     bz[nb]=$4;
     }
   }
if ( nconfig >= stepmin && nconfig <=stepmax )
  {
   print nconfig
   nstep=nstep+1;
#Find the Oxygen atom chemically bonded to each Hydrogen atom.
#This is effctive only for liquid water.
  if ( bonded == 1 && atoma == "H" && atomb == "O" )
    {
  print "CONFIGURE ",nconfig >FILENAME".HB"
  for (i=1;i<=na;i++)
    {
    RHP[i]=100.0;
     for (j=1;j<=nb;j++)
       {
       RHP0=(ax[i]-bx[j])^2+(ay[i]-by[j])^2+(az[i]-bz[j])^2
       if (RHP0 < RHP[i])
         {
         HPx[i]=bx[j];
         HPy[i]=by[j];
         HPz[i]=bz[j];
         RHP[i]=RHP0
         } 
       } 
    }}
# The main program to count the distance.
   for (ia=1;ia<=na;ia++)
       {
   for (ix=-1;ix<=1;ix++)
       {
   for (iy=-1;iy<=1;iy++)
       {
   for (iz=-1;iz<=1;iz++)
       {
# The normal count of pair distribution function.
for (ib=1;ib<=nb;ib++)
{
  RAB=sqrt((ax[ia]-bx[ib]-cellx*ix)^2+(ay[ia]-by[ib]-celly*iy)^2+(az[ia]-bz[ib]-cellz*iz)^2) ;
  RPH=sqrt((ax[ia]-HPx[ia])^2+(ay[ia]-HPy[ia])^2+(az[ia]-HPz[ia])^2) ;
  RPP=sqrt((HPx[ia]-bx[ib]-cellx*ix)^2+(HPy[ia]-by[ib]-celly*iy)^2+(HPz[ia]-bz[ib]-cellz*iz)^2) ;
  Rbin=int(RAB/deltar+0.5)
# The bonded pair distribution function.
  if (( bonded == 1 ) && ( RAB <= 5.0 ))
  {
   xmid=(ax[ia]+bx[ib]+cellx*ix)/2.0;
   ymid=(ay[ia]+by[ib]+celly*iy)/2.0;
   zmid=(az[ia]+bz[ib]+cellz*iz)/2.0;
   RR=0.99/4.0*RAB*RAB;
   nincore=0;
   for (iixx=1;iixx<=3;iixx++){ixx=ncell[iixx];if (nincore > 0) break;
   for (iiyy=1;iiyy<=3;iiyy++){iyy=ncell[iiyy];if (nincore > 0) break;
   for (iizz=1;iizz<=3;iizz++){izz=ncell[iizz];if (nincore > 0) break;
   for (iatom=1;iatom<=natom;iatom++){if (nincore > 0) break;
       RRcenter=(xmid-atomx[iatom]-ixx*cellx)^2+(ymid-atomy[iatom]-iyy*celly)^2+(zmid-atomz[iatom]-izz*cellz)^2;
       if ( RRcenter < RR ) nincore=nincore+1;
                              }}}} 
  }
  if ( Rbin <= ntotal ) pair00[Rbin]=pair00[Rbin]+1;
  if ( (Rbin <= ntotal) && (nincore == 0) && ( RAB <=5.0) ) {
    pair11[Rbin]=pair11[Rbin]+1;
    if ( atoma == "H" && atomb == "O" ){
    HBcos=(RAB*RAB+RPH*RPH-RPP*RPP)/2/RAB/RPH
    HBsin=sqrt(1-HBcos*HBcos+0.0000001)
    HBangle=atan2(HBsin,HBcos)*180.0/3.1415826
    if ( HBangle > 10.0 ){
    printf "%3s%3s%3s%3s%3s%3s%3s%12.6f%12.6f%12.6f\n",\
    ia,atoma,ib,atomb,ix,iy,iz,RAB,RPP,HBangle > FILENAME".HB"
    }}}}
# End of nb cycles
       }}}}
# end of na cycles
    }
# end of steps
    }}

END{
ncoord=0;
ncoord1=0;
Rho0=atomb/cellx/celly/cellz;
print "The below is the pair correlation function between ",atoma," and ", atomb,"." > FILENAME"_"atoma"_"atomb".txt";
print "The cell size is x=",cellx,", y=",celly,"and z=",cellz,"." > FILENAME"_"atoma"_"atomb".txt";
print "The dynamic steps are from ",stepmin,"to",stepmax,"." > FILENAME"_"atoma"_"atomb".txt";
print "The inteval of distance (delta r, in Angstrom) is,"deltar",", "and the total step is ",ntotal,"." > FILENAME"_"atoma"_"atomb".txt";
print "If the above parameters have some thing wrong, please change your file 'pait.txt'to correct them."  > FILENAME"_"atoma"_"atomb".txt"; 
print "     r       g(r)      n(r)      g0(r)     n0(r)" > FILENAME"_"atoma"_"atomb".txt";
for (i=1;i<=ntotal;i++)
  {
  ncoord=ncoord+pair00[i];
  ncoord1=ncoord1+pair11[i];
  printf  "%8.3f%10.3f%10.3f%10.3f%10.3f\n",i*deltar,\
    pair00[i]/(na*nb/cellx/celly/cellz)/(4*3.1415926*(i*deltar)^2*deltar)/nstep,ncoord/na/nstep, \
    pair11[i]/(na*nb/cellx/celly/cellz)/(4*3.1415926*(i*deltar)^2*deltar)/nstep,ncoord1/na/nstep \
 > FILENAME"_"atoma"_"atomb".txt";
}}
