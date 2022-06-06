/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "function/Function.h"
#include "function/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/File.h"

#include <cmath>

using namespace std;

namespace PLMD{
namespace function{

//+PLUMEDOC FUNCTION PATHCV
/*
This function is a path collective variable defined onto a set of other collective variables 
(rather than RMSDs to reference structures). 

This variable computes the progress ("s" component) along a parameterized curve 
(or "string of nodes") placed in the space of N other collective variables provided 
in the input. Similar as with the PATHMSD, a distance from the path can be 
defined (the "z" component), which is here the length of an N-dimensional vector. 
It can be used to explore the neighborhood of a path using a bias, but alternatively 
it can be used internally to optimize the path towards the mean transition path
(see \cite ensing12 ).


\par Examples
The following input tells plumed to print two torsion values and the progress
along a path defined in this 2-dimensional torsion space. The path is read from
the file "path.input". The format of the path input file is equal to that of the
path output, which can be generated as shown in the next example.
\verbatim
t1: TORSION ATOMS=5,7,9,15 
t2: TORSION ATOMS=7,9,15,17 
PATHCV LABEL=pcv ARG=t1,t2 INFILE=path.input
PRINT ARG=t1,t2,pcv.s STRIDE=10  FILE=colvar
\endverbatim
(See also \ref PRINT and \ref TORSION).

To generate an initial path as a straight interpolation between two free energy minima 
(defined as positions in the space of the collective variable arguments) use the GENPATH 
keyword. The first 3 numbers are the number of trailing nodes at the start of the path,
the number of nodes between the start and end, and the number of trailing nodes at the end.
The following example defines a straight path of 20 nodes between TORSIONs positions 
(-1.44,1.27) and (1.23,-1.21). At both sides, extra 10 trailing nodes are supplied. The
path is stored in output file "path.out", which can be used to create an input file.  
\verbatim
t1: TORSION ATOMS=5,7,9,15 
t2: TORSION ATOMS=7,9,15,17 
PATHCV LABEL=pcv ARG=t1,t2 GENPATH=10,20,10,-1.44,1.27,1.23,-1.21 FIXED=11,30 OUTFILE=path.out STRIDE=1
\endverbatim

The path can be optimized towards the highest transition density (or flux). The average crossing 
distance from the path is accumulated for each node. Every PACE md steps, the nodes which have
been passed (i.e. some average distance has been measured) are moved to reduce this average
distance. In the following example, a MOVINGRESTRAINT is used to pull the system along the path.
In addition, a RESTRAINT on the z-component keeps the sampling close to the path (aka the "tube"
picture). The "tube" potential may help converging the path, but note that it will affect a free
energy calculation along the path.
\verbatim
t1: TORSION ATOMS=5,7,9,15 
t2: TORSION ATOMS=7,9,15,17 
PATHCV LABEL=pcv ARG=t1,t2 INFILE=path.input PACE=100 OUTFILE=path.out
MOVINGRESTRAINT ARG=pcv.s LABEL=restraint STEP0=0 AT0=0.0 KAPPA0=1000.0 STEP1=10000 AT1=1.0 KAPPA1=1000.0
RESTRAINT ARG=pcv.z LABEL=tube KAPPA=50.0 AT=0.0
PRINT ARG=t1,t2,pcv.s STRIDE=10  FILE=colvar
\endverbatim
(See also \ref MOVINGRESTRAINT and \ref RESTRAINT).

The path is stored every STRIDE md steps in OUTFILE, in blocks separated by two empty lines
for easy reading by gnuplot. The first line of each block contains the time step. Use for example
the following gnuplot commands to visualize, for a 2-dimensional CV space, the sampled points from the
colvar file (see \ref PRINT) together with the optimized path:
\verbatim
> gnuplot
p 'colvar' every 10 u 2:3 w p lw 1 lt 2 pt 6
datafile = 'path.out'
stats datafile
rep for [IDX=1:STATS_blocks:10] datafile index (IDX-1) u 2:3 w lp t columnheader(1)
\endverbatim

The path can also be optimized by running metadynamics on the s-component. We refer to this scheme as
"path-metadynamics" (PMD). The implementation is robust, as the metadynamics Gaussians overwrite and
converge the free energy profile while the path is optimized.
\verbatim
t1: TORSION ATOMS=5,7,9,15 
t2: TORSION ATOMS=7,9,15,17 
PATHCV LABEL=pcv ARG=t1,t2 INFILE=path.input PACE=100 OUTFILE=path.out
METAD LABEL=metadyn ARG=pcv.s SIGMA=0.10 HEIGHT=0.2 PACE=100
RESTRAINT ARG=pcv.z LABEL=tube KAPPA=50.0 AT=0.0
PRINT ARG=t1,t2,pcv.s STRIDE=10  FILE=colvar
\endverbatim
HALFLIFE

To run a multiple-walker PMD simulation, we use the same syntax as in the file-based multiple-walker version of metadynamics.
Currently, this code only supports file-based communication for the path (that can be combined either with the file-based
or MPI-based version of multiple-walkers metadynamics). In a later version, we will implement MPI communication for the path.
\verbatim
# change WALKER_ID accordingly for each plumed.$WALKER_ID.dat file 
PATHCV LABEL=pcv ARG=t1,t2 INFILE=path.input OUTFILE=PATH PACE=100 STRIDE=100 WALKERS_RSTRIDE=100 WALKERS_ID=0 WALKERS_N=8 WALKERS_DIR=.
METAD LABEL=metadyn ARG=pcv.s SIGMA=0.10 HEIGHT=0.2 PACE=100 WALKERS_MPI
\endverbatim
In this example, the individual contributions of each walker to the path are stored in PATH.$WALKER_ID files. 
The actual shared path is printed in PATH_.$WALKER_ID files as is the same for all walkers.

It can be convenient to bias the sampling not only in the direction along the path, but also perpendicular to it. 
In the following example, both the s- and z-components of a fixed path are used as collective variables 
in a metadynamics calculation to compute the free energy along the path and to search for other transition paths.
\verbatim
t1: TORSION ATOMS=5,7,9,15 
t2: TORSION ATOMS=7,9,15,17 
PATHCV LABEL=pcv ARG=t1,t2 INFILE=pathcv.input
METAD LABEL=metadyn ARG=pcv.s,pcv.z SIGMA=0.10,0.10 HEIGHT=0.2 PACE=100
PRINT ARG=t1,t2,pcv.s,pcv.z STRIDE=10  FILE=colvar
\endverbatim
(See also \ref METAD and \ref PRINT).

*/
//+ENDPLUMEDOC

#define TOL 0.000001

class PathCV :  public Function{

private:
  unsigned ncv, nnodes;
  int pace, halflife, stride, istep;
  std::vector<int> fixednodes;
  std::vector< std::vector<double> > path, disp;
  std::vector<double> wsum, _wsum, scale;
  std::vector<double> genpath;
  double fadefact;
  std::string inputfile, outputfile;
  unsigned mw_n_;
  string mw_dir_;
  unsigned mw_id_;
  int mw_rstride_;
  vector<IFile*> ifiles;
  vector<string> ifilesnames;
  vector<string> ifilesnames_;
  OFile pathOfile;
  OFile pathOfile_;

  bool  readPath(IFile *ifile, vector<vector<double> > &path,vector<vector<double> > &disp,vector<double> &wsum,int &istep);
  void  generatePath();
  void  updatePath();
  void  reparamPath();
  void  writePath();
  void  writePath_();
  void  readMultipleWalkers();

public:
  PathCV(const ActionOptions&);
  ~PathCV();
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(PathCV,"PATHCV")

void PathCV::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG"); keys.use("PERIODIC");
  keys.add("optional","INFILE","an input file from which an initial path is read");
  keys.add("optional","GENPATH","to generate an initial path as a linear interpolation between 2 fixed nodes, supply 3 integers for the number of head, middle, and tail nodes, followed by the positions of the two fixed nodes");
  keys.add("optional","FIXED","the numbers of the two path nodes that remain fixed during a path optimization. If omitted, by default the first and last nodes are fixed");
  //keys.add("optional","SCALE","list of scaling factors for all CVs to stretch or shrink the CV space in certain directions in the case that different types of CVs (with different ranges) are used, such as a distance that runs from 0-10 Angstrom and a coordination number that runs from 0-1. Reasonable numbers are typically obtined by taking one divided by the (components of) the distance between the two fixed nodes. By default the scaling factors are set to 1.0");
  keys.add("optional","PACE","the path update pace (default=0)");
  keys.add("optional","OUTFILE","an output file in which the path updates are stored (default=PATH)");
  keys.add("optional","STRIDE","the frequency of writing the path to output file (default: equal to PACE)");
  keys.add("optional","HALFLIFE","the number of MD steps after which a previously measured path distance weighs only 50% in the average. This option may increase convergence by allowing to \"forget\" the memory of a bad initial guess path. A negative number sets it to infinity (=default)");
  keys.add("optional","WALKERS_ID", "walker id");
  keys.add("optional","WALKERS_N", "number of walkers");
  keys.add("optional","WALKERS_DIR", "shared directory with the path files from all the walkers");
  keys.add("optional","WALKERS_RSTRIDE","stride for reading hills files");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("s","default","the position on the path");
  keys.addOutputComponent("z","default","the distance from the path");
}


PathCV::~PathCV(){
  pathOfile.close();
  pathOfile_.close();
  // close files
  for(int i=0;i<mw_n_;++i){
   if(ifiles[i]->isOpen()) ifiles[i]->close();
   delete ifiles[i];
  }
}

PathCV::PathCV(const ActionOptions&ao):
Action(ao),
Function(ao),
nnodes(0),
pace(0),
halflife(-1),
stride(0),
outputfile("PATH"),
// Multiple walkers initialization
mw_n_(1), mw_dir_("./"), mw_id_(0), mw_rstride_(1)
{
  // open path output file
  parse("OUTFILE",outputfile);
  //  ofile.link(*this); 
  //  ofile.open(outputfile);

  // Multiple walkers
  parse("WALKERS_N",mw_n_);
  parse("WALKERS_ID",mw_id_);
  if(mw_n_<=mw_id_) error("walker ID should be a numerical value less than the total number of walkers");
  parse("WALKERS_DIR",mw_dir_);
  parse("WALKERS_RSTRIDE",mw_rstride_);
  for(int i=0;i<mw_n_;++i){
    string fname;
    string fname_;
    if(mw_n_>1) {
      stringstream out; out << i;
      fname = mw_dir_+"/"+outputfile+"."+out.str();
      fname_ = mw_dir_+"/"+outputfile+"_."+out.str();
    } else {
      fname = outputfile;
    }
    IFile *ifile = new IFile();
    ifile->link(*this);
    ifiles.push_back(ifile);                                                             
    ifilesnames.push_back(fname);
    if(mw_n_>1) ifilesnames_.push_back(fname_);
    if(ifile->FileExist(fname)){
      ifile->open(fname);
      //   if(plumed.getRestart()){
      //      log.printf("  Restarting from %s:",ifilesnames[i].c_str());                  
      //      readPath(ifiles[i]);                                                  
      //    }
      ifiles[i]->reset(false);
      // close only the walker own hills file for later writing
      if(i==mw_id_) ifiles[i]->close();
    }
  }

  //  open path output file for writing
  pathOfile.link(*this);
  pathOfile.enforceSuffix(""); 
  pathOfile.open(ifilesnames[mw_id_]);
  pathOfile.setHeavyFlush(); 

  if(mw_n_>1){
  pathOfile_.link(*this);
  pathOfile_.enforceSuffix(""); 
  pathOfile_.open(ifilesnames_[mw_id_]);
  pathOfile_.setHeavyFlush(); 
  }

  //  open the path input file
  parse("INFILE",inputfile);
  if(inputfile.size()!=0){
    IFile *ifp = new IFile();
    ifp->open(inputfile.c_str());
    int istep;
    if(readPath(ifp,path,disp,wsum,istep)){
     log.printf("  reading initial path from file: step=%d\n",istep);
    }else{
      log.printf("  At step=%d : Failed to read path file\n");
    }
    ifp->close();
  }


  // generate node positions in a straight line if no path was read from file
  parseVector("GENPATH",genpath);
  if( (genpath.size()==0) && (inputfile.size()==0) ){
    error("Please supply GENPATH or INFILE to initialize a path");
  }else if(genpath.size()!=0){
    if(inputfile.size()!=0){
      error("Found conflicting GENPATH and INFILE to initialize the path; please choose one");
    }
    generatePath();

    // allocate arrays to store average distance to path and weights and initialize with zero
    disp.resize(nnodes);
    for(unsigned i=0;i<nnodes;++i) disp[i].resize(ncv);
    wsum.resize(nnodes);
    for(unsigned i=0;i<nnodes;++i){
      for(unsigned j=0;j<ncv;++j) disp[i][j]=0.0;
      wsum[i] = 0.0;
    }
  }


  // fix two nodes when supplied in input, otherwise fix first and last node. Check if supplied nodes exist.
  parseVector("FIXED",fixednodes);
  if(fixednodes.size()==0){
    fixednodes.resize(2);
    fixednodes[0]=1;
    fixednodes[1]=nnodes;
  }
  if(fixednodes.size()!=2)
    error("The FIXED keyword requires two integers indicating the fixed nodes");
  fixednodes[0]--;
  fixednodes[1]--;
  if((fixednodes[0]<0)||(fixednodes[0]>=(int)nnodes)) error("first fixed node does not exist");
  if((fixednodes[1]<0)||(fixednodes[1]>=(int)nnodes)) error("second fixed node does not exist");

  // get the scaling factors for the CVs, or set to 1.0 if not supplied
  //parseVector("SCALE",scale);
  if(scale.size()==0){
    scale.resize(ncv);
    for(unsigned i=0;i<ncv;++i) scale[i]=1.0;
  }else if(scale.size()!=ncv){
    error("Size of SCALE array should be the same as the number of CVs as arguments");
  }

  // reparametrize the path (could be necessary if the path was read from file) 
  reparamPath();

  parse("PACE",pace);
  if(pace<0 ) error("frequency for path updates cannot be negative");
  stride=pace;   // Default: write path to file after each path update
  parse("HALFLIFE",halflife);
  if(halflife<0){
    fadefact=1.0;
  }else{
    //    fadefact = exp(log(0.5)/((double)halflife));  /* ah, "log" is reserved for something else...  */
    fadefact = exp(-0.693147180559945 / ( (double)halflife ));
  }
  parse("STRIDE",stride);
  if(stride<0 ) error("frequency for writing the path output cannot be negative");


  // print some feedback to log file
  log.printf("  path nodes:\n");
  for(unsigned inode=0;inode<path.size();inode++){
    for(unsigned icv=0;icv<path[0].size();icv++) log.printf("    %f",path[inode][icv]);
    log.printf("\n");
  }
  if(pace==0){
    log.printf("  path remains fixed\n");
  }else{
    log.printf("  path is flexible\n");
    log.printf("  path update pace: %d\n",pace);
    log.printf("  fixed path nodes:\n");
    log.printf("  scaling factors: %f",scale[0]);
    for(unsigned icv=1;icv<ncv;icv++) log.printf(", %f",scale[icv]);
    log.printf("\n");
    for(unsigned ii=0;ii<2;ii++){
      log.printf("    node [%d]: ",fixednodes[ii]+1);
      for(unsigned icv=0;icv<path[0].size();icv++) log.printf(" %f",path[fixednodes[ii]][icv]);
      log.printf("\n");
    }
    log.printf("  halflife of data: %d  (fade factor= %g)\n",halflife,fadefact);
  }
  if(stride>0){
    log.printf("  path is written to output file: %s\n",outputfile.c_str());
    log.printf("  frequency for writing: %d\n",stride);
    if(stride<pace) log.printf("WARNING: writing of path to file is more frequent than path is updated\n");
  }
  if(mw_n_>1){
    log.printf("  multiple walkers active: %d\n",mw_n_);
    log.printf("  walker id: %d\n",mw_id_);
    log.printf("  reading stride: %d\n",mw_rstride_);
    log.printf("  directory with path files: %s\n",mw_dir_.c_str());
  }

  addComponentWithDerivatives("s"); componentIsNotPeriodic("s");
  addComponentWithDerivatives("z"); componentIsNotPeriodic("z");
  checkRead();
  log<<"  Bibliography "<<plumed.cite("Díaz Leines and Ensing, Phys. Rev. Lett. 109, 020601 (2012)")<<"\n";
  log<<"  Bibliography "<<plumed.cite("Pérez de Alba Ortíz, Tiwari, Puthenkalathil and Ensing, J. Chem. Phys. 149, 072320 (2018)")<<"\n";
}


  bool PathCV::readPath(IFile *ifile, vector<vector<double> > &path,vector<vector<double> > &disp,vector<double> &wsum,int &istep){
    int idum,nn;
  double dummy;
  // read and check header line
   //log.printf("  readPath called\n");
   //if((ifile->FieldExist("Step")) && (ifile->FieldExist("Time")) && (ifile->FieldExist("nCV")) && (ifile->FieldExist("nNodes"))){
   if(ifile->FieldExist("Step")) {
   //log.printf("  step exists\n");
   if(ifile->FieldExist("Time")); //log.printf(" time exists\n");
   else return false;
   if(ifile->FieldExist("nCV")); //log.printf("  ncv exists\n"); 
   else return false;
   if(ifile->FieldExist("nNodes")); //log.printf("  nnodes exists\n");
   else return false;

    //log.printf("  sleeping for a sec\n");
    //sleep(1);
    //log.printf("  let's scan\n");
    ifile->scanField("Step",istep);
    //log.printf("  Read header step\n");
    ifile->scanField("Time",dummy);
    //log.printf("  Read header time\n");
    ifile->scanField("nCV",dummy);  ncv=(int)dummy;
    //log.printf("  Read header ncv\n");
    ifile->scanField("nNodes",dummy); nn=(int)dummy;
    //log.printf("  Read header nnodes\n");
    ifile->scanField();
    if( ncv != static_cast<unsigned>(getNumberOfArguments()) )
      error("Dimension of path in inputfile header should be equal to the number of arguments");
    if( nnodes==0 ){
      nnodes=nn;
    }else if( nn != (int)nnodes){ 
      error("Number of nodes has changed (multiple walkers with different paths lengths?)");
    }
    if(path.size() != nnodes){
      // allocate arrays to store path, average distance to path, and weights
      path.resize(nnodes);
      for(unsigned inode=0;inode<nnodes;++inode) path[inode].resize(ncv);
      disp.resize(nnodes);
      for(unsigned inode=0;inode<nnodes;++inode) disp[inode].resize(ncv);
      wsum.resize(nnodes);
    }
    // read the path nodes
    for(unsigned inode=0;inode<nnodes;++inode){
      if(ifile->FieldExist("node")) ifile->scanField("node",idum);
      for(unsigned icv=0;icv<ncv;++icv){
	if(!ifile->FieldExist(getPntrToArgument(icv)->getName()))
	  error("Unknown CV names in FIELDS header of path input file. Please fix.");
	//path[inode][icv]=0.0;
	ifile->scanField(getPntrToArgument(icv)->getName(),path[inode][icv]);
      }
      for(unsigned icv=0;icv<ncv;++icv){
	if(!ifile->FieldExist("z_avg_"+getPntrToArgument(icv)->getName()))
	  error("Unknown (z_avg_) CV names in FIELDS header of path input file. Please fix.");
	//disp[inode][icv]=0.0;
	ifile->scanField("z_avg_"+getPntrToArgument(icv)->getName(),disp[inode][icv]);
      }
      //wsum[inode]=0.0;
      if(ifile->FieldExist("wsum"))ifile->scanField("wsum",wsum[inode]);
      else return false;
      ifile->scanField();
    }
    // read last 3 lines
    std::string line;
    ifile->getline(line);
    ifile->getline(line);
    ifile->getline(line);
    return true;
  }else{ 
    return false;
  }
}


void PathCV::generatePath(){
  ncv = static_cast<unsigned>(getNumberOfArguments());
  if(genpath.size()!=2*ncv+3){
    char line[128];
    sprintf(line,"GENPATH requires 3 integers (numbers of head, middle and tail nodes) plus 2x %d values for the fixed node positions %ld",ncv,inputfile.size());
    error(line);
  }
  int nheadnodes=(int)genpath[0];
  int nmiddlenodes=(int)genpath[1];
  int ntailnodes=(int)genpath[2];
  nnodes = nheadnodes + nmiddlenodes + ntailnodes;
  path.resize(nnodes);
  for(unsigned inode=0;inode<nnodes;++inode) path[inode].resize(ncv);
  for(unsigned icv=0;icv<ncv;++icv){
    double dz =  (genpath[3+ncv+icv]-genpath[3+icv]) / (double) (nmiddlenodes-1);
    for(int inode=0;inode<nnodes;inode++){
      log.printf("  1: %f\n",genpath[3+icv]);
      log.printf("  2: %f\n",(double)(inode-nheadnodes));
      log.printf("  3: %f\n",dz);
      path[inode][icv] = genpath[3+icv] + (double)(inode-nheadnodes) * dz;
      log.printf("  node: %f\n",path[inode][icv]);
    }
  }
}


void PathCV::calculate(){
  int inodemin1 = 0, inodemin2 = 1;
  double distmin1 = 99999., distmin2 = 99999., dist;
  std::vector<double> v1(ncv);
  std::vector<double> v2(ncv);
  std::vector<double> v3(ncv);
  std::vector<double> p(ncv);
  for(unsigned icv=0;icv<ncv;++icv) p[icv]=getArgument(icv);

  // determine closest and second closest path nodes to current position
  for(unsigned inode=0;inode<nnodes;++inode){
    dist = 0.;
    for(unsigned icv=0;icv<ncv;++icv){
      dist += pow(scale[icv]*(getArgument(icv) - path[inode][icv]),2);
    }
    if( dist < distmin1 ){
      distmin2 = distmin1;
      inodemin2 = inodemin1;
      distmin1 = dist;
      inodemin1 = inode;
    }else if( dist < distmin2){
      distmin2 = dist;
      inodemin2 = inode;      
    }
  }
  int isign = inodemin1 - inodemin2;
  if(isign>1){
    isign=1;
  }else if(isign<-1){
    isign=-1;
  }
  inodemin2 = inodemin1 - isign;
  int inodemin3 = inodemin1 + isign;

  // compute some vectors with respect to closest nodes
  for(unsigned icv=0;icv<ncv;++icv){
    v1[icv] = scale[icv]*(path[inodemin1][icv] - p[icv]);
    v3[icv] = scale[icv]*(p[icv] - path[inodemin2][icv]);
  }
  if( (inodemin3 < 0) || (inodemin3 >= (int)nnodes)){
    for(unsigned icv=0;icv<ncv;++icv){
      v2[icv] = scale[icv]*(path[inodemin1][icv] - path[inodemin2][icv]);
    }
  }else{
    for(unsigned icv=0;icv<ncv;++icv){
      v2[icv] = scale[icv]*(path[inodemin3][icv] - path[inodemin1][icv]);
    }
  }

  // and some dot products
  double v1v1, v2v2, v3v3, v1v2;
  v1v1 = v2v2 = v3v3 = v1v2 = 0.;
  for(unsigned icv=0;icv<ncv;++icv){
    v1v1 += v1[icv] * v1[icv];
    v1v2 += v1[icv] * v2[icv];
    v2v2 += v2[icv] * v2[icv];
    v3v3 += v3[icv] * v3[icv];
  }

  // compute path_s value
  int nmiddle;
  double root, dx, path_s;
  Value* val_s_path=getPntrToComponent("s");
  root = sqrt( v1v2*v1v2 - v2v2 * ( v1v1 - v3v3) ); 
  dx = 0.5 * ( (root - v1v2) / v2v2 - 1.);
  path_s = inodemin1 + isign * dx - fixednodes[0];
  nmiddle = fixednodes[1] - fixednodes[0];
  path_s /= (double)nmiddle; 
  val_s_path->set(path_s);

  // compute derivatives of path_s
  double fact = isign * 0.5 / (v2v2 * (double)nmiddle );
  for(unsigned icv=0;icv<ncv;++icv){
    setDerivative(val_s_path, icv, fact/scale[icv]*( v2[icv] + (v2v2 * (v1[icv] + v3[icv]) - v1v2*v2[icv] )/root ) );
  }

  // compute vector v1 of current position to path and accumulate 
  // for the two closest nodes the (weighted) distance   
  double weight2 = -1.*dx;
  double weight1 = 1.0 + dx;
  double path_z = 0.0;
  if(weight1 > 1.0){                      /* current position is in front of first node */
    weight1 = 1.0; weight2 = 0.0;
  }else if(weight2 > 1.0){                /* current position is after last node */
    weight1 = 0.0; weight2 = 1.0;
  }
  for(unsigned icv=0;icv<ncv;++icv){
    v1[icv] = p[icv] - path[inodemin1][icv] - dx * (path[inodemin1][icv] - path[inodemin2][icv]);
    disp[inodemin1][icv] += weight1 * v1[icv];
    disp[inodemin2][icv] += weight2 * v1[icv];
    path_z += v1[icv] * v1[icv];
  }
  wsum[inodemin1] *= fadefact;
  wsum[inodemin2] *= fadefact;
  wsum[inodemin1] += weight1;
  wsum[inodemin2] += weight2;

  // set value and derivatives of path_z
  Value* val_z_path=getPntrToComponent("z");
  path_z=sqrt(path_z); val_z_path->set(path_z);
  for(unsigned icv=0;icv<ncv;++icv){
    setDerivative(val_z_path, icv, v1[icv]/path_z );
  }

  // update path nodes after pace steps and not on step 0
  if( (pace>0) && (getStep()>0) && (getStep()%pace==0) ){
    updatePath();  
    reparamPath();
  }

  // write path to output file
  if( (stride>0) && (getStep()%stride==0) ){
    writePath();
  }

  // read paths from multiple walkers after mw_rstride steps and not on step 0
  if( (mw_n_>1) && (getStep()>0) && (getStep()%mw_rstride_==0) ){
    readMultipleWalkers();
    reparamPath();
  }

  // write path to output file
  if( (mw_n_>1) && (stride>0) && (getStep()%stride==0) ){
    writePath_();
  }

}


void PathCV::updatePath(){

  // set weights of fixed nodes to zero 
  int inow=getStep();
  log.printf("  At step=%d : Weight of fixed node 0 of walker %d is %f\n",inow,mw_id_,wsum[fixednodes[0]]);
  log.printf("  At step=%d : Weight of fixed node 1 of walker %d is %f\n",inow,mw_id_,wsum[fixednodes[1]]);
  wsum[fixednodes[0]] = 0;
  wsum[fixednodes[1]] = 0;

  // update node if wsum>0 (I could consider to update only if wsum>TOL??)
  for(unsigned inode=0;inode<nnodes;++inode){
    if(wsum[inode]>0){
      for(unsigned icv=0;icv<ncv;++icv){
	path[inode][icv] += disp[inode][icv] / wsum[inode];
	disp[inode][icv] = 0.0;
      }
    }
  }
}


/* ------------------------------------------------------------------------------------ */
/* Simple algorithm to make the nodes that define the path more equidistant.            */
/* Adapted from J. Chem. Phys. 125, 024106 (2006) String method, Vanden Eijnden et al.  */
/*  - first the nodes between the fixed nodes are redistributed, then the end pieces    */
/*    are done, such that their separation is equal to that found for the middle piece. */
void PathCV::reparamPath(){

  double dx, dr, lenpiece;
  std::vector<double> len(nnodes), sumlen(nnodes), sfrac(nnodes);
  std::vector< std::vector<double> > newpath;
  newpath.resize(nnodes);
  for(unsigned i=0;i<nnodes;++i) newpath[i].resize(ncv);

  int istart = fixednodes[0];
  int iend = fixednodes[1];
  int nmiddle = iend - istart + 1;

  /* compute actual distances between path nodes */
  len[istart] = sumlen[istart] = 0.0;
  for(int i=istart+1;i<iend+1;i++){
    dr = 0.; 
    for(unsigned j=0;j<ncv;j++){
      dx = scale[j]*(path[i][j]-path[i-1][j]);
      dr += dx*dx;
    }
    len[i] = sqrt(dr);
    sumlen[i] = sumlen[i-1] + len[i];
  }


  /* iterate until some tolerance */
  unsigned iter = 0;
  double prevsum = 0.;
  while( fabs(sumlen[iend] - prevsum) > TOL){
    prevsum = sumlen[iend];

    /* compute cumulative target distances between path nodes */
    dr = sumlen[iend];
    for(int i=0;i<nmiddle;i++){
      sfrac[istart+i] = (double)i*dr/(double)(nmiddle-1);
    }

    /* compute new nodes positions */
    for(int i=istart+1;i<iend;i++){
      int k = istart;
      while( !((sumlen[k] < sfrac[i]) && (sumlen[k+1] >= sfrac[i])) ){
        k++;
        if(k >= iend+1){
          error("Error in reparametrizing path");
        }
      }
      dr = (sfrac[i]-sumlen[k])/len[k+1];
      for(unsigned j=0;j<ncv;j++){  
        newpath[i][j] = path[k][j] + dr*(path[k+1][j]-path[k][j]) ;
      }
    }

    /* copy new path */
    for(int i=istart+1;i<iend;i++){
      for(unsigned j=0;j<ncv;j++) path[i][j] = newpath[i][j];
    }

    /* compute actual distances between path nodes */
    len[istart] = sumlen[istart] = 0.0;
    for(int i=istart+1;i<iend+1;i++){
      dr = 0.; 
      for(unsigned j=0;j<ncv;j++){
	dx = scale[j]*(path[i][j]-path[i-1][j]);
	dr += dx*dx;
      }
      len[i] = sqrt(dr);
      sumlen[i] = sumlen[i-1] + len[i];
    }

    iter++;
  }
  lenpiece = sumlen[iend]/(nmiddle-1);

  /* now same thing for the path ends */
  double tfrac = sfrac[iend]/(nmiddle-1);
  istart = fixednodes[1];
  iend = nnodes - 1;
  nmiddle = iend - istart + 1;

  /* compute actual distances between path nodes */
  len[istart] = sumlen[istart] = 0.;
  for(int i=istart+1;i<iend+1;i++){
    dr = 0.; 
    for(unsigned j=0;j<ncv;j++){
      dx = scale[j]*(path[i][j]-path[i-1][j]);
      dr += dx*dx;
    }
    len[i] = sqrt(dr);
    sumlen[i] = sumlen[i-1] + len[i];
    sfrac[i] = (double)(i-istart) * tfrac; 
  }

  /* iterate until some tolerance */
  unsigned iter2 = 0;
  prevsum = 0.;
  while( fabs(sumlen[iend] - prevsum) > TOL){
    prevsum = sumlen[iend];
    
    /* compute new points */
    for(int i=istart+1;i<iend+1;i++){
      int k = istart;
      while( !((sumlen[k] < sfrac[i]) && (sumlen[k+1] >= sfrac[i])) ){
        k++;
        if(k >= iend){
	  //          error("Error in reparametrizing path");
	  k = iend-1;
	  break;
        }
      }
      dr = (sfrac[i]-sumlen[k])/len[k+1];
      for(unsigned j=0;j<ncv;j++){  
        newpath[i][j] = path[k][j] + dr*(path[k+1][j]-path[k][j]) ;
      }
    }
  
    /* copy new points */
    for(int i=istart+1;i<iend+1;i++){
      for(unsigned j=0;j<ncv;j++) path[i][j] = newpath[i][j];
    }

    /* compute actual distances between path nodes */
    len[istart] = sumlen[istart] = 0.;
    for(int i=istart+1;i<iend+1;i++){
      dr = 0.; 
      for(unsigned j=0;j<ncv;j++){
        dx = scale[j]*(path[i][j]-path[i-1][j]);
        dr += dx*dx;
      }
      len[i] = sqrt(dr);
      sumlen[i] = sumlen[i-1] + len[i];
    }

    iter2++;
  }


  /* now same thing for the path other end */
  istart = fixednodes[0];
  iend = 0;
  nmiddle = istart - iend + 1;
  
  /* compute actual distances between path nodes */
  len[istart] = sumlen[istart] = 0.;
  for(int i=istart-1;i>iend-1;i--){
    dr = 0.; 
    for(unsigned j=0;j<ncv;j++){
      dx = scale[j]*(path[i][j]-path[i+1][j]);
      dr += dx*dx;
    }
    len[i] = sqrt(dr);
    sumlen[i] = sumlen[i+1] + len[i];
    sfrac[i] = (double)(istart-i) * tfrac; 
  }

  /* iterate until some tolerance */
  unsigned iter1 = 0;
  prevsum = 0.;
  while( fabs(sumlen[iend] - prevsum) > TOL){
    prevsum = sumlen[iend];
    
    /* compute new points */
    for(int i=istart-1;i>iend-1;i--){
      int k = istart;
      while( !((sumlen[k] < sfrac[i]) && (sumlen[k-1] >= sfrac[i])) ){
        k--;
        if(k <= iend){
//          error("Error in reparametrizing path");
            k = iend+1;
            break;
            
        }
      }
      dr = (sfrac[i]-sumlen[k])/len[k-1];
      for(unsigned j=0;j<ncv;j++){  
        newpath[i][j] = path[k][j] + dr*(path[k-1][j]-path[k][j]) ;
      }
    }
  
    /* copy new points */
    for(int i=istart-1;i>iend-1;i--){
      for(unsigned j=0;j<ncv;j++) path[i][j] = newpath[i][j];
    }
  
    /* compute actual distances between path nodes */
    len[istart] = sumlen[istart] = 0.;
    for(int i=istart-1;i>iend-1;i--){
      dr = 0.; 
      for(unsigned j=0;j<ncv;j++){
        dx = scale[j]*(path[i][j]-path[i+1][j]);
        dr += dx*dx;
      }
      len[i] = sqrt(dr);
      sumlen[i] = sumlen[i+1] + len[i];
    }

    iter1++;
  }
  log.printf("  path nodes redistributed in %d + %d + %d cycles. Length=%lf\n",iter1,iter,iter2,lenpiece);

}


void PathCV::writePath(){
  pathOfile.fmtField(" %f");
  pathOfile.printField("Step",(int)getStep());
  pathOfile.printField("Time",getTime());
  pathOfile.printField("nCV",(int)ncv);  
  pathOfile.printField("nNodes",(int)nnodes);  
  pathOfile.printField();
  for(unsigned inode=0;inode<nnodes;inode++){
    pathOfile.printField("node",(int)inode); // first colum, number of node
    for(unsigned icv=0;icv<ncv;icv++){
      pathOfile.printField(getPntrToArgument(icv)->getName(),path[inode][icv]); 
    }
    for(unsigned icv=0;icv<ncv;icv++){
      pathOfile.printField("z_avg_"+getPntrToArgument(icv)->getName(),disp[inode][icv]); 
    }
    pathOfile.printField("wsum",wsum[inode]); 
    pathOfile.printField(); 
  }
  pathOfile.printField();
  pathOfile.printField();
  pathOfile.flush();
}

void PathCV::writePath_(){
  pathOfile_.fmtField(" %f");
  pathOfile_.printField("Step",(int)getStep());
  pathOfile_.printField("Time",getTime());
  pathOfile_.printField("nCV",(int)ncv);  
  pathOfile_.printField("nNodes",(int)nnodes);  
  pathOfile_.printField();
  for(unsigned inode=0;inode<nnodes;inode++){
    pathOfile_.printField("node",(int)inode); // first colum, number of node
    for(unsigned icv=0;icv<ncv;icv++){
      pathOfile_.printField(getPntrToArgument(icv)->getName(),path[inode][icv]); 
    }
    for(unsigned icv=0;icv<ncv;icv++){
      pathOfile_.printField("z_avg_"+getPntrToArgument(icv)->getName(),disp[inode][icv]); 
    }
    pathOfile_.printField("wsum",wsum[inode]); 
    pathOfile_.printField(); 
  }
  pathOfile_.printField();
  pathOfile_.printField();
  pathOfile_.flush();
}


void PathCV::readMultipleWalkers(){

  // declaring variables
  std::vector< std::vector<double> > mw_path, mw_disp, tot_path,_tot_path;
  std::vector<double> mw_wsum;
  int inow=getStep();
  bool read_;

  if(tot_path.size() != nnodes){
    tot_path.resize(nnodes);
    _tot_path.resize(nnodes);

    for(unsigned inode=0;inode<nnodes;++inode) tot_path[inode].resize(ncv);
    for(unsigned inode=0;inode<nnodes;++inode) _tot_path[inode].resize(ncv);

  }

  /* print for debugging
  log.printf("  current path nodes:\n");
  for(unsigned inode=0;inode<path.size();inode++){
    for(unsigned icv=0;icv<path[0].size();icv++) log.printf("    %f",path[inode][icv]);
    log.printf("    %f",wsum[inode]);
    log.printf("\n");
  } */

  // first, multiply current path nodes by current weights
  log.printf("  I've done my change to the path\n");
  for(unsigned inode=0;inode<nnodes;++inode){
    for(unsigned icv=0;icv<ncv;++icv) {
      tot_path[inode][icv] = path[inode][icv] * wsum[inode];
      _tot_path[inode][icv] = path[inode][icv];
      //log.printf("    %f",path[inode][icv]);
    }
    //log.printf("    %f",wsum[inode]);
    //log.printf("\n"); 
  }
  // initialize array for sum of weights of all walkers
  _wsum.resize(nnodes);
  for(unsigned inode=0;inode<nnodes;++inode){
   _wsum[inode] = wsum[inode];
  }

  // next, read path from each walker and add it to the total multipied by its weights
  for(int i=0;i<mw_n_;++i){
    if(i==mw_id_) continue;
    // if the file is not open yet 
    if(!(ifiles[i]->isOpen())){
      // check if it exists now and open it!
      if(ifiles[i]->FileExist(ifilesnames[i])) {
	ifiles[i]->open(ifilesnames[i]);
	ifiles[i]->reset(false);
      }else{
	log.printf("  At step=%d : Walker %d has not produced path output yet\n",inow,i);
	continue;
      }
    }

    // otherwise read the path
    // barriers to sync walkers
    comm.Barrier();
    multi_sim_comm.Barrier();
    // first read to skip initial guess paths
    if(inow == mw_rstride_) {
     readPath(ifiles[i],mw_path,mw_disp,mw_wsum,istep); 
     ifiles[i]->reset(false);
    }
    // actual reading
    read_=false;
    while(!read_){
      log.printf("  At step=%d : Attempting to read path from walker %d\n",inow,i);
      read_=readPath(ifiles[i],mw_path,mw_disp,mw_wsum,istep);
      //log.printf("  From file %s\n",ifilesnames[i].c_str());
      //log.printf("  Address: %p\n",(void *)ifiles[i]);
      //std::string line;
      //ifiles[i]->getline(line);
      //log.printf("  Line: %s\n",line.c_str());
      ifiles[i]->reset(false);
    }
      log.printf("  At step=%d : Reading path[step=%d] from walker %d\n",inow,istep,i);
      for(unsigned inode=0;inode<nnodes;++inode){
	for(unsigned icv=0;icv<ncv;++icv){
	  tot_path[inode][icv] += mw_wsum[inode]*mw_path[inode][icv];
          _tot_path[inode][icv] += mw_path[inode][icv];
          //log.printf("    %f",mw_path[inode][icv]);
	}
	_wsum[inode] += mw_wsum[inode];
        //log.printf("    %f",mw_wsum[inode]);
        //log.printf("\n");
      }

  }
  // finally, divide again the path nodes by the total weights and reparametrize nodes
  for(unsigned inode=0;inode<nnodes;++inode){
    if(_wsum[inode] > 0) for(unsigned icv=0;icv<ncv;++icv) path[inode][icv] = tot_path[inode][icv] / _wsum[inode];
    else for(unsigned icv=0;icv<ncv;++icv) path[inode][icv] = _tot_path[inode][icv] / mw_n_;
  }

  ///* print for debugging */ /*
   //log.printf("  new path nodes:\n");
   //for(unsigned inode=0;inode<path.size();inode++){
     //for(unsigned icv=0;icv<path[0].size();icv++) log.printf("    %f",path[inode][icv]);
     //log.printf("    %f",wsum[inode]);
     //log.printf("\n");
   //} */

  reparamPath();

  ///* print for debugging */ /*readpa
   //log.printf("  new path nodes reparametrized:\n");
   //for(unsigned inode=0;inode<path.size();inode++){
     //for(unsigned icv=0;icv<path[0].size();icv++) log.printf("    %f",path[inode][icv]);
     //log.printf("    %f",wsum[inode]);
     //log.printf("\n");
   //} */

}


}
}


