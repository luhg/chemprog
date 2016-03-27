##############################################################
BEGIN{ 
  simcutoff=0.98;
  ngood=0;
  issimilar=0;
  nfirst=0;
  nstruct=0;
  }
##################
{
if ( nfirst == 0 ) {natom=$1;nfirst=1;}
if ( $1 == natom ){
nstruct=nstruct+1;
getline;
mycomment=$0;
for (i=1;i<=natom; i++){
  getline;
  nsym[i]=$1;
  atomx[i]=$2;
  atomy[i]=$3;
  atomz[i]=$4;
}
  usrindex();
  if ( ngood > 1){
    issimilar=0;
    for (ib=1;ib<=ngood;ib++){
      msim=0.0;
      for (ic=1;ic<=12;ic++){
        msim=msim+abs(usr[ic]-musr[ib,ic]);
        }
      similar=1/(1+msim/12);
      # print similar,simcutoff;
      if ( similar > simcutoff ) { issimilar=1;break;}
    }
  }
  if ( issimilar < 1 ) {    
    ngood=ngood+1;   
    printf "%-10s\n",natom;
    print mycomment;
    for (j=1;j<=natom;j++)
       {
       printf "%-10s %10.3f %10.3f %10.3f\n",nsym[j],atomx[j],atomy[j],atomz[j];
       }
    en[ngood]=mycomment;
    for (j=1;j<=12;j++){musr[ngood,j]=usr[j];}
    }else {print similar,simcutoff,nstruct,ngood,mycomment,en[ngood];}   
}}

##########################################################
function usrindex(){
# ctd centroid
cx=0.0;
cy=0.0;
cz=0.0;
for (ia=1;ia<=natom;ia++){
  cx=cx+atomx[ia];
  cy=cy+atomy[ia];
  cz=cz+atomz[ia];
  }
ctdx=cx/natom;
ctdy=cy/natom;
ctdz=cz/natom;
# cst,fct atom closest to ctd (cst), atom farthest to ctd (fct),
cstr0=100;
fctr0=0;
for (ia=1;ia<=natom;ia++) {
  usrr=dist(atomx[ia],atomy[ia],atomz[ia],ctdx,ctdy,ctdz);
  if (usrr < cstr0) {
     cstr0=usrr;
     ncst=ia
  }
  if (usrr > fctr0) {
     fctr0=usrr;
     nfct=ia;
  }
}
cstx=atomx[ncst];
csty=atomy[ncst];
cstz=atomz[ncst];
fctx=atomx[nfct];
fcty=atomy[nfct];
fctz=atomz[nfct];
# ftf atom farthest to fct (ftf).
ftfr0=0;
for (ia=1;ia<=natom;ia++) {
  usrr=dist(atomx[ia],atomy[ia],atomz[ia],fctx,fcty,fctz);
  if (usrr > ftfr0) {
     ftfr0=usrr;
     nftf=ia
  }
}
ftfx=atomx[nftf];
ftfy=atomy[nftf];
ftfz=atomz[nftf];
# distance
for (ia=1;ia<=natom;ia++) {
  ctdr[ia]=dist(atomx[ia],atomy[ia],atomz[ia],ctdx,ctdy,ctdz);
  cstr[ia]=dist(atomx[ia],atomy[ia],atomz[ia],cstx,csty,cstz);
  fctr[ia]=dist(atomx[ia],atomy[ia],atomz[ia],fctx,fcty,fctz);
  ftfr[ia]=dist(atomx[ia],atomy[ia],atomz[ia],ftfx,ftfy,ftfz);
}
#means
ctdrall=0.0;cstrall=0.0;fctrall=0.0;ftfrall=0.0;
for (ia=1;ia<=natom;ia++) {
  ctdrall=ctdrall+ctdr[ia];
  cstrall=cstrall+cstr[ia];
  fctrall=fctrall+fctr[ia];
  ftfrall=ftfrall+ftfr[ia];
}
usr[1]=ctdrall/natom;
usr[2]=cstrall/natom;
usr[3]=fctrall/natom;
usr[4]=ftfrall/natom;
ctd2=0.0;ctd3=0.0;
cst2=0.0;cst3=0.0;
fct2=0.0;fct3=0.0;
ftf2=0.0;ftf3=0.0;
for (ia=1;ia<=natom;ia++) {
  ctd2=ctd2+(ctdr[ia]-usr[1])^2;
  ctd3=ctd3+(ctdr[ia]-usr[1])^3;
  cst2=cst2+(cstr[ia]-usr[2])^2;
  cst3=cst3+(cstr[ia]-usr[2])^3;
  fct2=fct2+(fctr[ia]-usr[3])^2;
  fct3=fct3+(fctr[ia]-usr[3])^3;
  ftf2=ftf2+(ftfr[ia]-usr[4])^2;
  ftf3=ftf3+(ftfr[ia]-usr[4])^3;
}
usr[5]=sqrt(ctd2/natom);
usr[6]=sqrt(cst2/natom);
usr[7]=sqrt(fct2/natom);
usr[8]=sqrt(ftf2/natom);
usr[9] =abs(ctd3/natom)^(1/3)*sign(ctd3/natom);
usr[10]=abs(cst3/natom)^(1/3)*sign(cst3/natom);
usr[11]=abs(fct3/natom)^(1/3)*sign(fct3/natom);
usr[12]=abs(ftf3/natom)^(1/3)*sign(ftf3/natom);
}

function dist(x1,y1,z1,x2,y2,z2){ 
mydist=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
return mydist;
}

function abs(x) {
  return x < 0 ? -x : x
}

function sign(x) {
  return x < 0 ? -1 : 1
}
