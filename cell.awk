#（1）读入abc或ABC，分数坐标或直角坐标
#（2）转换为abc和ABC，分数坐标和直角坐标
#（3）根据需要输出指定格式的abc或ABC，分数坐标或直角坐标
# 变量说明：
# natom：原子总数；ntype[]:各种原子总数；
# input，output：输入和输出格式；
# ABC：晶胞矢量(ax,ay,az,bx,by,bz,cx,cy,cz);
# abc：晶胞参数(aa,bb,cc,alpha,beta,gamma)；
# 直接读入的直接坐标：atomcx,atomcy,atomcz,atomsymbol;
# 直接读入的分数坐标：atomfx,atomfy,atomfz,atomsymbol;
# 重新排序后直角坐标：atcx,atcy,atcz,atsymbol;
# 重新排序后分数坐标：atfx,atfy,atfz,atsymbol;

BEGIN{
PI=3.1415926;
natom=0;
symbol="Aa Bb Cc Dd Ee Ff Gg Hh Ii Jj Kk Ll Mm Nn Oo Pp Qq Rr Ss Tt Uu Vv Ww Xx Yy Zz"
split(symbol,mysymbol);
if (input == "" ) {input="pdb";}
if (output == "" ) {output="cif";}
}
# Read input file
(input == "pdb" ){Labc=1;Lcart=1;readpdb();}
(input == "gen" ){LABC=1;Lfract=1;readgen();}
(input == "fract" ) {Labc=1;Lfract=1;readfract();}
(input == "cif" ) {Labc=1;Lfract=1;readcif();}
(input == "poscar" ) {LABC=1;Lfract=1;readposcar();}

END{
if (LABC == 1) {ABC2abc();}
if (Labc == 1) {abc2ABC();}
if (Lfract == 1) {fract2cart();}
if (Lcart == 1) {cart2fract();}
reorderatom();
if (output == "poscar") {writeposcar();}
if (output == "cif") {writecif();}
if (output == "findsym") {writefindsym();}
}

# Read functions: readpdb, readcif, readfract, readposcar
function readpdb(){
for (i=1;i<=1000;i++) {
  getline;
  if ($1 == "CRYST1" ){
    aa=$2;
    bb=$3;
    cc=$4;
    alpha=$5;
    beta=$6;
    gamma=$7;
    break;
    }
  }
# Read Catersian coordinates
  for (i=1;i<=10000;i++){
    getline;
    if ($1 != "TER" ) {
      natom=natom+1;
      atomsymbol[natom]=$11;
      atomcx[natom]=$6;
      atomcy[natom]=$7;
      atomcz[natom]=$8;
      }else{
        exit;}
      }
  exit;
}

function readfract(){
# Read a b c alpha beta gamma
getline;
aa=$1;
bb=$2;
cc=$3;
alpha=$4;
beta=$5;
gamma=$6;
# Read fractional coordinates
for (i=1;i<=10000;i++){
  getline;
  if (NF == 4 ) {
    natom=natom+1;
    atomsymbol[natom]=$1;
    atomfx[natom]=$2;
    atomfy[natom]=$3;
    atomfz[natom]=$4;
    }else{
      exit;}
    }
exit;
}

function readcif(){
for (i=1;i<=1000;i++) {
# Read a b c alpha beta gamma
  getline;
  if ($1 == "_cell_length_a" ){ aa=$2;}
  if ($1 == "_cell_length_b" ){ bb=$2;}
  if ($1 == "_cell_length_c" ){ cc=$2;}
  if ($1 == "_cell_angle_alpha" ){ alpha=$2;}
  if ($1 == "_cell_angle_beta" ){ beta=$2;}
  if ($1 == "_cell_angle_gamma" ){ gamma=$2;}
# Read fractional coordinates
  if ($1 == "_atom_site_occupancy" ){ 
    nr0=0;
    for (j=1;j<=10000;j++){
      getline;
      if (NR == nr0) {break;}else {nr0=NR}
      if (NF == 8 ) {
        natom=natom+1;
        atomsymbol[natom]=$2;
        atomfx[natom]=$3;
        atomfy[natom]=$4;
        atomfz[natom]=$5;
        }else{
          exit;}
        }
      exit;
     }
  }
}

function readposcar(){
# Read scale, A, B, C
getline;
cellscale=$1;
getline;
ax=$1*cellscale;
ay=$2*cellscale;
az=$3*cellscale;
getline;
bx=$1*cellscale;
by=$2*cellscale;
bz=$3*cellscale;
getline;
cx=$1*cellscale;
cy=$2*cellscale;
cz=$3*cellscale;
# Read element symbols and numbers
getline;
if ( "ABCDEFGHIJKLMNOPQRSTUVWXYZ" ~ toupper(substr($1,1,1))  ) { 
  for (i=1;i<=NF;i++) {atomtype[i]=$i;}
  getline;
  }
for (i=1;i<=NF;i++) {
  if (atomtype[i] == "" ) {
    atomtype[i]=mysymbol[i];
  }
  for (j=1;j<=$i;j++) {
    natom=natom+1;
    atomsymbol[natom]=atomtype[i];
  }
}
getline;
# Select Dynamics
if (toupper(substr($1,1,1)) == "S") {getline;}
# Direct coordinates
if (toupper(substr($1,1,1)) == "D") {
  Lfract=1;
  for (i=1;i<=natom;i++) {
    getline;
    atomfx[i]=$1;
    atomfy[i]=$2;
    atomfz[i]=$3;
  }
}
# Cartesian coordinates
else {
  Lcart=1;
  for (i=1;i<=natom;i++) {
    getline;
    atomcx[i]=$1;
    atomcy[i]=$2;
    atomcz[i]=$3;
  }
}
exit;
}
#read Gen of DFTB
function readgen(){
  natom=$1;
  myfract=$2;
  getline;
  for (i=1;i<=NF;i++){
    nsymbol[i]=$i;
    }
  for (i=1;i<=natom;i++){
    getline;
    atomsymbol[i]=nsymbol[$2];
    atomfx[i]=$3;
    atomfy[i]=$4;
    atomfz[i]=$5;
    }
  getline;
  getline;ax=$1;ay=$2;az=$3;
  getline;bx=$1;by=$2;bz=$3;
  getline;cx=$1;cy=$2;cz=$3;
  exit;
}

# Convert functions: abc2ABC, ABC2abc, fract2cart, cart2frac, reorderatom

function abc2ABC(){ 
# Convert the a b c to vectors A B C.
  ax=aa;
  ay=0;
  az=0;
  alpha0=alpha/180*PI;
  beta0=beta/180*PI;
  gamma0=gamma/180*PI;
  bx=bb*cos(gamma0);
  by=bb*sin(gamma0);
  bz=0;
  cx=cc*cos(beta0);
  cy=cc*(cos(alpha0)-cos(beta0)*cos(gamma0))/sin(gamma0);
  cz=cc*sqrt(1+2*cos(alpha0)*cos(beta0)*cos(gamma0)-cos(alpha0)^2-cos(beta0)^2-cos(gamma0)^2)/sin(gamma0);
}

function ABC2abc() {
# Convert the vectors A B C to a b c and alpha beta gamma.
  aa=sqrt(ax^2+ay^2+az^2);
  bb=sqrt(bx^2+by^2+bz^2);
  cc=sqrt(cx^2+cy^2+cz^2);
  ab=sqrt((ax-bx)^2+(ay-by)^2+(az-bz)^2);
  ac=sqrt((ax-cx)^2+(ay-cy)^2+(az-cz)^2);
  bc=sqrt((bx-cx)^2+(by-cy)^2+(bz-cz)^2);
  cosa=-(bc^2-bb^2-cc^2)/2/bb/cc;
  cosb=-(ac^2-aa^2-cc^2)/2/aa/cc;
  cosc=-(ab^2-aa^2-bb^2)/2/aa/bb;
  alpha=atan2(sqrt(1-cosa^2),cosa)/PI*180.0;
  beta=atan2(sqrt(1-cosb^2),cosb)/PI*180.0;
  gamma=atan2(sqrt(1-cosc^2),cosc)/PI*180.0;
}

function fract2cart(i){ 
  for (i=1;i<=natom;i++) {
    atomcx[i]=ax*atomfx[i];
    atomcy[i]=bx*atomfx[i]+by*atomfy[i];
    atomcz[i]=cx*atomfx[i]+cy*atomfy[i]+cz*atomfz[i];
  }
}

function cart2fract(){
# Convert matrix from Catersian to fractional coordinates 
  vol=ax*by*cz;
  acx=by*cz/vol;
  acy=-bx*cz/vol;
  acz=(bx*cy-by*cx)/vol;
  bcx=0;
  bcy=ax*cz/vol;
  bcz=-ax*cy/vol;
  ccx=0;
  ccy=0;
  ccz=ax*by/vol;
  for (i=1;i<=natom;i++) {
    atomfx[i]=acx*atomcx[i];
    atomfy[i]=bcx*atomcx[i]+bcy*atomcy[i];
    atomfz[i]=ccx*atomcx[i]+ccy*atomcy[i]+ccz*atomcz[i];
  }
}

function reorderatom(i){ 
nat=0;
# Reorder all atoms in cell
for (i=1;i<=natom;i++) {
  ntype[atomsymbol[i]]=ntype[atomsymbol[i]]+1;
}
for (i in ntype )
  {
  for (j=1; j<=natom; j++)
    {
    if (atomsymbol[j]==i) {
      nat=nat+1;
      atsymbol[nat]=i;
      atfx[nat]=atomfx[j];
      atfy[nat]=atomfy[j];
      atfz[nat]=atomfz[j];
      atcx[nat]=atomcx[j];
      atcy[nat]=atomcy[j];
      atcz[nat]=atomcz[j];
       }
    }
  } 
}


# Write functions
function writeposcar(i) {
  printf "paw_pbe" 
  for (i in ntype ) {
  printf "%4s",i;
  }
  print "\n1.000" 
  printf "%10.5f%10.5f%10.5f\n",ax,ay,az;
  printf "%10.5f%10.5f%10.5f\n",bx,by,bz;
  printf "%10.5f%10.5f%10.5f\n",cx,cy,cz;
  for (i in ntype ) {
  printf "%4s",i;
  }
  printf "\n";
  for (i in ntype) {
  printf "%4s",ntype[i];
  }
  printf "\n";
  print "Direct";
  for (i=1;i<=natom;i++) {
  printf "%10.5f%10.5f%10.5f   !%2s\n",atfx[i],atfy[i],atfz[i],atsymbol[i];
  }
}

function writecif(i){
print "data_cif" ;
print "_audit_creation_date              2012-12-1" ;
print "_audit_creation_method            'CIF Generate by PMMP'";
print "_symmetry_space_group_name_H-M    'P1'";
print "_symmetry_Int_Tables_number       1";
print "_symmetry_cell_setting            triclinic";
print "loop_";
print "_symmetry_equiv_pos_as_xyz";
print "  x,y,z";
printf "_cell_length_a                  %10.4f\n",aa;
printf "_cell_length_b                  %10.4f\n",bb;
printf "_cell_length_c                  %10.4f\n",cc;
printf "_cell_angle_alpha               %10.4f\n",alpha;
printf "_cell_angle_beta                %10.4f\n",beta;
printf "_cell_angle_gamma               %10.4f\n",gamma;
print "loop_";
print "_atom_site_label";
print "_atom_site_type_symbol";
print "_atom_site_fract_x";
print "_atom_site_fract_y";
print "_atom_site_fract_z";
print "_atom_site_U_iso_or_equiv";
print "_atom_site_adp_type";
print "_atom_site_occupancy";
for (i=1;i<=natom;i++){
  col1=atsymbol[i]i;
  printf "%-7s%-4s%10.5f%10.5f%10.5f  0.00000   Uiso    1.00\n",col1,atsymbol[i],atfx[i],atfy[i],atfz[i];
}
}

function writefindsym(i) {
  print "input for findsym\n 0.1\n 1";
  printf "%10.5f%10.5f%10.5f\n",ax,ay,az;
  printf "%10.5f%10.5f%10.5f\n",bx,by,bz;
  printf "%10.5f%10.5f%10.5f\n",cx,cy,cz;
  print " 1\n 1 0 0 \n 0 1 0 \n 0 0 1";
  print natom;
  for (i in ntype ) {
    for ( j=1;j<=ntype[i];j++){
      printf "%3s",i;
    }
  }
  printf "\n"; 
  for (i=1;i<=natom;i++) {
  printf "%10.5f%10.5f%10.5f\n",atfx[i],atfy[i],atfz[i];
  }
}

