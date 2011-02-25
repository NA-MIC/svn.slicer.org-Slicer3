/*
 *  sh2hot.cxx
 *  
 *
 *  Created by Antonio Tristán Vega on 17/09/10.
 *  Copyright 2010 E.T.S.I de Telecomunicación. All rights reserved.
 *
 *  The code in this file is used to translate between spherical harmonics and
 *  higher order tensors. Note that this is only intended (the conversion is only
 *  possible) with HOT of even order.
 */
 


#ifndef _sh2hot_cxx
#define _sh2hot_cxx

#include "sphericalHarmonics.h"
// These two header files are used to solve a system of nonlinear, polynomial
// equations of three equations (corresponding to each component of the 
// Z-eigenvector) with three unknowns (corresponding to each coordinate x, y, z):
#include <vnl/algo/vnl_rpoly_roots.h>
#include <vnl/vnl_real_npolynomial.h>
#include <vnl/algo/vnl_rnpoly_solve.h>

namespace shmaths
{
   /**
    *  This function implements the precomputed transition matrices between the coefficients of a real,
    *  symmetric function defined over the unit sphere in the basis of Spherical Harmonics and the
    *  coefficients of this same function in the basis of higher order tensors. These two basis are
    *  known to be basis of the same functional space. The details on this implementation may be 
    *  found in:
    *
    *      Maxime Descoteaux, Elaine Angelino, Shaun Fitzgibbons, and Rachid Deriche
    *      "Apparent Diffusion Coefficients from High Angular
    *                       Resolution Diffusion Images: Estimation and Applications"
    *      Tech. Report 5681, Institut National de Recherche en Informatique et Automatique
    *      September 2005
    *
    *  The matrix for basis change M as in eq. (38) is implemented for orders 0, 2, 4, 6, and 8.
    *  Orders above 8 are never considered in diffusion imaging, and hence this code should be
    *  useful for most of practical applications.
    *
    *  NOTE: The integrals in eq. (38) are NOT numerically evaluated; instead, analytical computation
    *  is performed with Maple 13 and then the result is numerically evaluated with arbitrary
    *  precission (15 digits). [See the file SH2HOT2.mw.]
    *
    *  Fortunatelly, the matrix M is highly sparse in all cases, so most of the values are 0.
    *
    */
   bool generateSH2HOTMatrix( const unsigned int order, SHMatrixType& M )
   {
      if(   order > 8   ||   order%2 != 0   )
         return false;
      M.set_size( (order+1)*(order+2)/2, (order+1)*(order+2)/2 );
      M.fill( 0.00000000000000000000f );
      switch( order ){
         case 0:
            M(0,0)=+3.544907701811032;
            break;
         case 2:
            M(0,0)=+1.181635900603677; M(0,2)=+1.181635900603677; M(0,5)=+1.181635900603677; 
            M(1,2)=-0.915291232863769; M(1,5)=+0.915291232863769; 
            M(2,3)=+1.830582465727538; 
            M(3,0)=+1.056887279361603; M(3,2)=-0.528443639680801; M(3,5)=-0.528443639680801; 
            M(4,1)=-1.830582465727538; 
            M(5,4)=+1.830582465727538;
            break;
         case 4:
            M(0,0)=+0.708981540362206; M(0,2)=+1.417963080724413; M(0,4)=+0.708981540362206; M(0,9)=+1.417963080724413; M(0,11)=+1.417963080724413; M(0,14)=+0.708981540362206; 
            M(1,2)=-0.784535342454659; M(1,4)=-0.784535342454659; M(1,9)=+0.784535342454659; M(1,14)=+0.784535342454659; 
            M(2,5)=+1.569070684909318; M(2,7)=+1.569070684909318; M(2,12)=+1.569070684909318; 
            M(3,0)=+0.905903382309945; M(3,2)=+0.452951691154973; M(3,4)=-0.452951691154973; M(3,9)=+0.452951691154973; M(3,11)=-0.905903382309945; M(3,14)=-0.452951691154973; 
            M(4,1)=-1.569070684909318; M(4,3)=-1.569070684909318; M(4,10)=-1.569070684909318; 
            M(5,6)=+1.569070684909318; M(5,8)=+1.569070684909318; M(5,13)=+1.569070684909318; 
            M(6,4)=+0.199732921787032; M(6,11)=-1.198397530722193; M(6,14)=+0.199732921787032; 
            M(7,7)=-1.694790041061752; M(7,12)=+0.564930013687251; 
            M(8,2)=-0.905903382309945; M(8,4)=+0.150983897051658; M(8,9)=+0.905903382309945; M(8,14)=-0.150983897051658; 
            M(9,5)=+0.854093899641589; M(9,7)=-0.640570424731192; M(9,12)=-0.640570424731192; 
            M(10,0)=+0.270088205852269; M(10,2)=-0.810264617556807; M(10,4)=+0.101283077194601; M(10,9)=-0.810264617556807; M(10,11)=+0.202566154389202; M(10,14)=+0.101283077194601; 
            M(11,1)=-0.854093899641589; M(11,3)=+0.640570424731192; M(11,10)=+0.640570424731192; 
            M(12,6)=+1.811806764619891; M(12,8)=-0.301967794103315; M(12,13)=-0.301967794103315; 
            M(13,3)=+0.564930013687251; M(13,10)=-1.694790041061752; 
            M(14,8)=-0.798931687148128; M(14,13)=+0.798931687148128;
            break;
         case 6:
            M(0,0)=+0.506415385973005; M(0,2)=+1.519246157919014; M(0,4)=+1.519246157919014; M(0,6)=+0.506415385973005; M(0,13)=+1.519246157919014; M(0,15)=+3.038492315838027; M(0,17)=+1.519246157919014; M(0,22)=+1.519246157919014; M(0,24)=+1.519246157919014; M(0,27)=+0.506415385973005; 
            M(1,2)=-0.653779452045549; M(1,4)=-1.307558904091099; M(1,6)=-0.653779452045549; M(1,13)=+0.653779452045549; M(1,17)=-0.653779452045549; M(1,22)=+1.307558904091099; M(1,24)=+0.653779452045549; M(1,27)=+0.653779452045549; 
            M(2,7)=+1.307558904091099; M(2,9)=+2.615117808182197; M(2,11)=+1.307558904091099; M(2,18)=+2.615117808182197; M(2,20)=+2.615117808182197; M(2,25)=+1.307558904091099; 
            M(3,0)=+0.754919485258288; M(3,2)=+1.132379227887432; M(3,6)=-0.377459742629144; M(3,13)=+1.132379227887432; M(3,17)=-1.132379227887432; M(3,24)=-1.132379227887432; M(3,27)=-0.377459742629144; 
            M(4,1)=-1.307558904091099; M(4,3)=-2.615117808182197; M(4,5)=-1.307558904091099; M(4,14)=-2.615117808182197; M(4,16)=-2.615117808182197; M(4,23)=-1.307558904091099; 
            M(5,8)=+1.307558904091099; M(5,10)=+2.615117808182197; M(5,12)=+1.307558904091099; M(5,19)=+2.615117808182197; M(5,21)=+2.615117808182197; M(5,26)=+1.307558904091099; 
            M(6,4)=+0.272363075164135; M(6,6)=+0.272363075164135; M(6,15)=-1.634178450984808; M(6,17)=-1.361815375820673; M(6,22)=+0.272363075164135; M(6,24)=-1.361815375820673; M(6,27)=+0.272363075164135; 
            M(7,9)=-2.311077328720572; M(7,11)=-2.311077328720572; M(7,18)=+0.770359109573524; M(7,20)=-1.540718219147048; M(7,25)=+0.770359109573524; 
            M(8,2)=-1.235322794059016; M(8,4)=-1.029435661715847; M(8,6)=+0.205887132343169; M(8,13)=+1.235322794059016; M(8,17)=+0.205887132343169; M(8,22)=+1.029435661715847; M(8,24)=-0.205887132343169; M(8,27)=-0.205887132343169; 
            M(9,7)=+1.164673499511258; M(9,9)=+0.291168374877814; M(9,11)=-0.873505124633443; M(9,18)=+0.291168374877814; M(9,20)=-1.747010249266887; M(9,25)=-0.873505124633443; 
            M(10,0)=+0.368302098889458; M(10,2)=-0.736604197778916; M(10,4)=-0.966793009584827; M(10,6)=+0.138113287083547; M(10,13)=-0.736604197778916; M(10,15)=-1.933586019169654; M(10,17)=+0.414339861250640; M(10,22)=-0.966793009584827; M(10,24)=+0.414339861250640; M(10,27)=+0.138113287083547; 
            M(11,1)=-1.164673499511258; M(11,3)=-0.291168374877814; M(11,5)=+0.873505124633443; M(11,14)=-0.291168374877814; M(11,16)=+1.747010249266887; M(11,23)=+0.873505124633443; 
            M(12,8)=+2.470645588118033; M(12,10)=+2.058871323431694; M(12,12)=-0.411774264686339; M(12,19)=+2.058871323431694; M(12,21)=-0.823548529372678; M(12,26)=-0.411774264686339; 
            M(13,3)=+0.770359109573524; M(13,5)=+0.770359109573524; M(13,14)=-2.311077328720572; M(13,16)=-1.540718219147048; M(13,23)=-2.311077328720572; 
            M(14,10)=-1.089452300656539; M(14,12)=-1.089452300656539; M(14,19)=+1.089452300656539; M(14,26)=+1.089452300656539; 
            M(15,6)=-0.045741696509789; M(15,17)=+0.686125447646828; M(15,24)=-0.686125447646828; M(15,27)=+0.045741696509789; 
            M(16,11)=+0.792269423793498; M(16,20)=-1.584538847586995; M(16,25)=+0.158453884758700; 
            M(17,4)=+0.337824817398841; M(17,6)=-0.033782481739884; M(17,15)=-2.026948904393045; M(17,17)=+0.168912408699420; M(17,22)=+0.337824817398841; M(17,24)=+0.168912408699420; M(17,27)=-0.033782481739884; 
            M(18,9)=-1.480274183795270; M(18,11)=+0.555102818923226; M(18,18)=+0.493424727931757; M(18,20)=+0.370068545948818; M(18,25)=-0.185034272974409; 
            M(19,2)=-0.493424727931757; M(19,4)=+0.493424727931757; M(19,6)=-0.030839045495735; M(19,13)=+0.493424727931757; M(19,17)=-0.030839045495735; M(19,22)=-0.493424727931757; M(19,24)=+0.030839045495735; M(19,27)=+0.030839045495735; 
            M(20,7)=+0.312069198822651; M(20,9)=-0.780172997056628; M(20,11)=+0.195043249264157; M(20,18)=-0.780172997056628; M(20,20)=+0.390086498528314; M(20,25)=+0.195043249264157; 
            M(21,0)=+0.068099082174715; M(21,2)=-0.510743116310365; M(21,4)=+0.383057337232774; M(21,6)=-0.021280963179599; M(21,13)=-0.510743116310365; M(21,15)=+0.766114674465548; M(21,17)=-0.063842889538796; M(21,22)=+0.383057337232774; M(21,24)=-0.063842889538796; M(21,27)=-0.021280963179599; 
            M(22,1)=-0.312069198822651; M(22,3)=+0.780172997056628; M(22,5)=-0.195043249264157; M(22,14)=+0.780172997056628; M(22,16)=-0.390086498528314; M(22,23)=-0.195043249264157; 
            M(23,8)=+0.986849455863514; M(23,10)=-0.986849455863514; M(23,12)=+0.061678090991470; M(23,19)=-0.986849455863514; M(23,21)=+0.123356181982939; M(23,26)=+0.061678090991470; 
            M(24,3)=+0.493424727931757; M(24,5)=-0.185034272974409; M(24,14)=-1.480274183795270; M(24,16)=+0.370068545948818; M(24,23)=+0.555102818923226; 
            M(25,10)=-1.351299269595363; M(25,12)=+0.135129926959536; M(25,19)=+1.351299269595363; M(25,26)=-0.135129926959536; 
            M(26,5)=-0.158453884758700; M(26,16)=+1.584538847586995; M(26,23)=-0.792269423793498; 
            M(27,12)=+0.274450179058731; M(27,21)=-0.914833930195771; M(27,26)=+0.274450179058731;
            break;
         case 8:
            M(0,0)=+0.393878633534559; M(0,2)=+1.575514534138236; M(0,4)=+2.363271801207355; M(0,6)=+1.575514534138236; M(0,8)=+0.393878633534559; M(0,17)=+1.575514534138236; M(0,19)=+4.726543602414709; M(0,21)=+4.726543602414709; M(0,23)=+1.575514534138236; M(0,30)=+2.363271801207355; M(0,32)=+4.726543602414709; M(0,34)=+2.363271801207355; M(0,39)=+1.575514534138236; M(0,41)=+1.575514534138236; M(0,44)=+0.393878633534559; 
            M(1,2)=-0.554721959311375; M(1,4)=-1.664165877934125; M(1,6)=-1.664165877934125; M(1,8)=-0.554721959311375; M(1,17)=+0.554721959311375; M(1,21)=-1.664165877934125; M(1,23)=-1.109443918622750; M(1,30)=+1.664165877934125; M(1,32)=+1.664165877934125; M(1,39)=+1.664165877934125; M(1,41)=+1.109443918622750; M(1,44)=+0.554721959311375; 
            M(2,9)=+1.109443918622750; M(2,11)=+3.328331755868251; M(2,13)=+3.328331755868251; M(2,15)=+1.109443918622750; M(2,24)=+3.328331755868251; M(2,26)=+6.656663511736502; M(2,28)=+3.328331755868251; M(2,35)=+3.328331755868251; M(2,37)=+3.328331755868251; M(2,42)=+1.109443918622750; 
            M(3,0)=+0.640537745067638; M(3,2)=+1.601344362669095; M(3,4)=+0.960806617601457; M(3,6)=-0.320268872533819; M(3,8)=-0.320268872533819; M(3,17)=+1.601344362669095; M(3,19)=+1.921613235202914; M(3,21)=-0.960806617601457; M(3,23)=-1.281075490135276; M(3,30)=+0.960806617601457; M(3,32)=-0.960806617601457; M(3,34)=-1.921613235202914; M(3,39)=-0.320268872533819; M(3,41)=-1.281075490135276; M(3,44)=-0.320268872533819; 
            M(4,1)=-1.109443918622750; M(4,3)=-3.328331755868251; M(4,5)=-3.328331755868251; M(4,7)=-1.109443918622750; M(4,18)=-3.328331755868251; M(4,20)=-6.656663511736502; M(4,22)=-3.328331755868251; M(4,31)=-3.328331755868251; M(4,33)=-3.328331755868251; M(4,40)=-1.109443918622750; 
            M(5,10)=+1.109443918622750; M(5,12)=+3.328331755868251; M(5,14)=+3.328331755868251; M(5,16)=+1.109443918622750; M(5,25)=+3.328331755868251; M(5,27)=+6.656663511736502; M(5,29)=+3.328331755868251; M(5,36)=+3.328331755868251; M(5,38)=+3.328331755868251; M(5,43)=+1.109443918622750; 
            M(6,4)=+0.293314080945991; M(6,6)=+0.586628161891982; M(6,8)=+0.293314080945991; M(6,19)=-1.759884485675947; M(6,21)=-2.933140809459912; M(6,23)=-1.173256323783965; M(6,30)=+0.293314080945991; M(6,32)=-2.933140809459912; M(6,34)=-2.933140809459912; M(6,39)=+0.586628161891982; M(6,41)=-1.173256323783965; M(6,44)=+0.293314080945991; 
            M(7,11)=-2.488852507852923; M(7,13)=-4.977705015705846; M(7,15)=-2.488852507852923; M(7,24)=+0.829617502617641; M(7,26)=-3.318470010470564; M(7,28)=-4.148087513088205; M(7,35)=+1.659235005235282; M(7,37)=-0.829617502617641; M(7,42)=+0.829617502617641; 
            M(8,2)=-1.330347624371248; M(8,4)=-2.438970644680622; M(8,6)=-0.886898416247499; M(8,8)=+0.221724604061875; M(8,17)=+1.330347624371248; M(8,21)=-0.886898416247499; M(8,23)=+0.443449208123749; M(8,30)=+2.438970644680622; M(8,32)=+0.886898416247499; M(8,39)=+0.886898416247499; M(8,41)=-0.443449208123749; M(8,44)=-0.221724604061875; 
            M(9,9)=+1.254263768704432; M(9,11)=+1.567829710880539; M(9,13)=-0.627131884352216; M(9,15)=-0.940697826528324; M(9,24)=+1.567829710880539; M(9,26)=-1.254263768704432; M(9,28)=-2.822093479584971; M(9,35)=-0.627131884352216; M(9,37)=-2.822093479584971; M(9,42)=-0.940697826528324; 
            M(10,0)=+0.396633029573262; M(10,2)=-0.396633029573262; M(10,4)=-1.834427761776338; M(10,6)=-0.892424316539840; M(10,8)=+0.148737386089973; M(10,17)=-0.396633029573262; M(10,19)=-3.668855523552677; M(10,21)=-2.677272949619521; M(10,23)=+0.594949544359893; M(10,30)=-1.834427761776338; M(10,32)=-2.677272949619521; M(10,34)=+0.892424316539840; M(10,39)=-0.892424316539840; M(10,41)=+0.594949544359893; M(10,44)=+0.148737386089973; 
            M(11,1)=-1.254263768704432; M(11,3)=-1.567829710880539; M(11,5)=+0.627131884352216; M(11,7)=+0.940697826528324; M(11,18)=-1.567829710880539; M(11,20)=+1.254263768704432; M(11,22)=+2.822093479584971; M(11,31)=+0.627131884352216; M(11,33)=+2.822093479584971; M(11,40)=+0.940697826528324; 
            M(12,10)=+2.660695248742497; M(12,12)=+4.877941289361244; M(12,14)=+1.773796832494998; M(12,16)=-0.443449208123749; M(12,25)=+4.877941289361244; M(12,27)=+3.547593664989996; M(12,29)=-1.330347624371248; M(12,36)=+1.773796832494998; M(12,38)=-1.330347624371248; M(12,43)=-0.443449208123749; 
            M(13,3)=+0.829617502617641; M(13,5)=+1.659235005235282; M(13,7)=+0.829617502617641; M(13,18)=-2.488852507852923; M(13,20)=-3.318470010470564; M(13,22)=-0.829617502617641; M(13,31)=-4.977705015705846; M(13,33)=-4.148087513088205; M(13,40)=-2.488852507852923; 
            M(14,12)=-1.173256323783965; M(14,14)=-2.346512647567929; M(14,16)=-1.173256323783965; M(14,25)=+1.173256323783965; M(14,29)=-1.173256323783965; M(14,36)=+2.346512647567929; M(14,38)=+1.173256323783965; M(14,43)=+1.173256323783965; 
            M(15,6)=-0.085384500151605; M(15,8)=-0.085384500151605; M(15,21)=+1.280767502274079; M(15,23)=+1.195383002122474; M(15,32)=-1.280767502274079; M(15,39)=+0.085384500151605; M(15,41)=-1.195383002122474; M(15,44)=+0.085384500151605; 
            M(16,13)=+1.478902924414529; M(16,15)=+1.478902924414529; M(16,26)=-2.957805848829057; M(16,28)=-1.478902924414529; M(16,35)=+0.295780584882906; M(16,37)=-2.662025263946152; M(16,42)=+0.295780584882906; 
            M(17,4)=+0.630606325811169; M(17,6)=+0.567545693230052; M(17,8)=-0.063060632581117; M(17,19)=-3.783637954867017; M(17,21)=-2.837728466150262; M(17,23)=+0.252242530324468; M(17,30)=+0.630606325811169; M(17,32)=-2.837728466150262; M(17,34)=+0.630606325811169; M(17,39)=+0.567545693230052; M(17,41)=+0.252242530324468; M(17,44)=-0.063060632581117; 
            M(18,11)=-2.763178476417838; M(18,13)=-1.726986547761149; M(18,15)=+1.036191928656689; M(18,24)=+0.921059492139279; M(18,26)=-1.151324365174099; M(18,28)=+1.726986547761149; M(18,35)=+0.575662182587050; M(18,37)=+0.345397309552230; M(18,42)=-0.345397309552230; 
            M(19,2)=-0.921059492139279; M(19,6)=+0.863493273880574; M(19,8)=-0.057566218258705; M(19,17)=+0.921059492139279; M(19,21)=+0.863493273880574; M(19,23)=-0.115132436517410; M(19,32)=-0.863493273880574; M(19,39)=-0.863493273880574; M(19,41)=+0.115132436517410; M(19,44)=+0.057566218258705; 
            M(20,9)=+0.582529171135615; M(20,11)=-0.873793756703423; M(20,13)=-1.092242195879279; M(20,15)=+0.364080731959760; M(20,24)=-0.873793756703423; M(20,26)=-2.184484391758557; M(20,28)=+1.092242195879279; M(20,35)=-1.092242195879279; M(20,37)=+1.092242195879279; M(20,42)=+0.364080731959760; 
            M(21,0)=+0.127118286726135; M(21,2)=-0.826268863719879; M(21,4)=-0.238346787611504; M(21,6)=+0.675315898232594; M(21,8)=-0.039724464601917; M(21,17)=-0.826268863719879; M(21,19)=-0.476693575223007; M(21,21)=+2.025947694697781; M(21,23)=-0.158897858407669; M(21,30)=-0.238346787611504; M(21,32)=+2.025947694697781; M(21,34)=-0.238346787611504; M(21,39)=+0.675315898232594; M(21,41)=-0.158897858407669; M(21,44)=-0.039724464601917; 
            M(22,1)=-0.582529171135615; M(22,3)=+0.873793756703423; M(22,5)=+1.092242195879279; M(22,7)=-0.364080731959760; M(22,18)=+0.873793756703423; M(22,20)=+2.184484391758557; M(22,22)=-1.092242195879279; M(22,31)=+1.092242195879279; M(22,33)=-1.092242195879279; M(22,40)=-0.364080731959760; 
            M(23,10)=+1.842118984278559; M(23,14)=-1.726986547761149; M(23,16)=+0.115132436517410; M(23,27)=-3.453973095522298; M(23,29)=+0.345397309552230; M(23,36)=-1.726986547761149; M(23,38)=+0.345397309552230; M(23,43)=+0.115132436517410; 
            M(24,3)=+0.921059492139279; M(24,5)=+0.575662182587050; M(24,7)=-0.345397309552230; M(24,18)=-2.763178476417838; M(24,20)=-1.151324365174099; M(24,22)=+0.345397309552230; M(24,31)=-1.726986547761149; M(24,33)=+1.726986547761149; M(24,40)=+1.036191928656689; 
            M(25,12)=-2.522425303244678; M(25,14)=-2.270182772920210; M(25,16)=+0.252242530324468; M(25,25)=+2.522425303244678; M(25,29)=+0.252242530324468; M(25,36)=+2.270182772920210; M(25,38)=-0.252242530324468; M(25,43)=-0.252242530324468; 
            M(26,5)=-0.295780584882906; M(26,7)=-0.295780584882906; M(26,20)=+2.957805848829057; M(26,22)=+2.662025263946152; M(26,31)=-1.478902924414529; M(26,33)=+1.478902924414529; M(26,40)=-1.478902924414529; 
            M(27,14)=+0.512307000909632; M(27,16)=+0.512307000909632; M(27,27)=-1.707690003032106; M(27,29)=-1.195383002122474; M(27,36)=+0.512307000909632; M(27,38)=-1.195383002122474; M(27,43)=+0.512307000909632; 
            M(28,8)=+0.010717813501466; M(28,23)=-0.300098778041036; M(28,34)=+0.750246945102590; M(28,41)=-0.300098778041036; M(28,44)=+0.010717813501466; 
            M(29,15)=-0.300098778041036; M(29,28)=+1.500493890205180; M(29,37)=-0.900296334123108; M(29,42)=+0.042871254005862; 
            M(30,6)=-0.109580580141874; M(30,8)=+0.007827184295848; M(30,21)=+1.643708702128115; M(30,23)=-0.109580580141874; M(30,32)=-1.643708702128115; M(30,39)=+0.109580580141874; M(30,41)=+0.109580580141874; M(30,44)=-0.007827184295848; 
            M(31,13)=+1.014519036400841; M(31,15)=-0.253629759100210; M(31,26)=-2.029038072801683; M(31,28)=+0.253629759100210; M(31,35)=+0.202903807280168; M(31,37)=+0.456533566380379; M(31,42)=-0.050725951820042; 
            M(32,4)=+0.281376954282889; M(32,6)=-0.168826172569733; M(32,8)=+0.007034423857072; M(32,19)=-1.688261725697331; M(32,21)=+0.844130862848666; M(32,23)=-0.028137695428289; M(32,30)=+0.281376954282889; M(32,32)=+0.844130862848666; M(32,34)=-0.070344238570722; M(32,39)=-0.168826172569733; M(32,41)=-0.028137695428289; M(32,44)=+0.007034423857072; 
            M(33,11)=-0.871814606355354; M(33,13)=+1.089768257944193; M(33,15)=-0.163465238691629; M(33,24)=+0.290604868785118; M(33,26)=+0.726512171962795; M(33,28)=-0.272442064486048; M(33,35)=-0.363256085981398; M(33,37)=-0.054488412897210; M(33,42)=+0.054488412897210; 
            M(34,2)=-0.214625919507681; M(34,4)=+0.536564798769203; M(34,6)=-0.201211799538451; M(34,8)=+0.006707059984615; M(34,17)=+0.214625919507681; M(34,21)=-0.201211799538451; M(34,23)=+0.013414119969230; M(34,30)=-0.536564798769203; M(34,32)=+0.201211799538451; M(34,39)=+0.201211799538451; M(34,41)=-0.013414119969230; M(34,44)=-0.006707059984615; 
            M(35,9)=+0.102610815720113; M(35,11)=-0.538706782530591; M(35,13)=+0.448922318775492; M(35,15)=-0.056115289846937; M(35,24)=-0.538706782530591; M(35,26)=+0.897844637550985; M(35,28)=-0.168345869540810; M(35,35)=+0.448922318775492; M(35,37)=-0.168345869540810; M(35,42)=-0.056115289846937; 
            M(36,0)=+0.017101802620019; M(36,2)=-0.239425236680263; M(36,4)=+0.448922318775492; M(36,6)=-0.149640772925164; M(36,8)=+0.004676274153911; M(36,17)=-0.239425236680263; M(36,19)=+0.897844637550985; M(36,21)=-0.448922318775492; M(36,23)=+0.018705096615646; M(36,30)=+0.448922318775492; M(36,32)=-0.448922318775492; M(36,34)=+0.028057644923468; M(36,39)=-0.149640772925164; M(36,41)=+0.018705096615646; M(36,44)=+0.004676274153911; 
            M(37,1)=-0.102610815720113; M(37,3)=+0.538706782530591; M(37,5)=-0.448922318775492; M(37,7)=+0.056115289846937; M(37,18)=+0.538706782530591; M(37,20)=-0.897844637550985; M(37,22)=+0.168345869540810; M(37,31)=-0.448922318775492; M(37,33)=+0.168345869540810; M(37,40)=+0.056115289846937; 
            M(38,10)=+0.429251839015363; M(38,12)=-1.073129597538406; M(38,14)=+0.402423599076902; M(38,16)=-0.013414119969230; M(38,25)=-1.073129597538406; M(38,27)=+0.804847198153805; M(38,29)=-0.040242359907690; M(38,36)=+0.402423599076902; M(38,38)=-0.040242359907690; M(38,43)=-0.013414119969230; 
            M(39,3)=+0.290604868785118; M(39,5)=-0.363256085981398; M(39,7)=+0.054488412897210; M(39,18)=-0.871814606355354; M(39,20)=+0.726512171962795; M(39,22)=-0.054488412897210; M(39,31)=+1.089768257944193; M(39,33)=-0.272442064486048; M(39,40)=-0.163465238691629; 
            M(40,12)=-1.125507817131554; M(40,14)=+0.675304690278933; M(40,16)=-0.028137695428289; M(40,25)=+1.125507817131554; M(40,29)=-0.028137695428289; M(40,36)=-0.675304690278933; M(40,38)=+0.028137695428289; M(40,43)=+0.028137695428289; 
            M(41,5)=-0.202903807280168; M(41,7)=+0.050725951820042; M(41,20)=+2.029038072801683; M(41,22)=-0.456533566380379; M(41,31)=-1.014519036400841; M(41,33)=-0.253629759100210; M(41,40)=+0.253629759100210; 
            M(42,14)=+0.657483480851246; M(42,16)=-0.046963105775089; M(42,27)=-2.191611602837486; M(42,29)=+0.109580580141874; M(42,36)=+0.657483480851246; M(42,38)=+0.109580580141874; M(42,43)=-0.046963105775089; 
            M(43,7)=+0.042871254005862; M(43,22)=-0.900296334123108; M(43,33)=+1.500493890205180; M(43,40)=-0.300098778041036; 
            M(44,16)=-0.085742508011725; M(44,29)=+0.600197556082072; M(44,38)=-0.600197556082072; M(44,43)=+0.085742508011725;
            break;
         default:
            break;
      }
      return true;
   }
   
   /** This function is intended to compute all possible combinations of powers os each component of the
    vector [x,y,z]^t to be included in each term of the tensor product. For a tensor of order L, the tensor
    sumation is defined as terms of the form:
    
             p)          q)         r)             p)            q)            r)
    T_{1,1, ..., 1,2,2, ...,2,3,3, ..., 3} · x·x· ... ·x · y·y· ... ·y · z·z· ... ·z = T_{i_1,i_2,...,i_L} · x^p · y^q · z^r
    
    where p+q+r=L. Since we deal only with hypersymmetric tensors, the term T_{i_1,i_2,...,i_L} is 
    invariant under any permutation of the set of indices {i_1,i_2,...,i_L}:
    
    T_{1,1,2,2,2,3,3} = T_{3,1,2,1,2,3,2} = T_{3,3,2,2,2,1,1} = ...
    
    It is easy to prove that there are exactly (L+1)(L+2)/2 independent terms for a symmetric tensor
    of order L. The purpose of this method is to find the powers of x, y, and z for each corresponding
    term in the summation.
    
    The memory to store the indices has to be allocated externally. Note that HOT form a basis for the
    same functional space as even-order SH, so the memory required can be allocated via:
    
    unsigned int* nx = new unsigned int[ getNumberOfEvenAssociatedLegendrePolynomials(order) ];
    */
   void computeSetOfHOTPowers( const unsigned int order, unsigned int* nx, unsigned int* ny, unsigned int* nz )
   {
      // Fill the buffers. The methodology is quite simple; since we need that nx+ny+nz=l for
      // each term, we may pick up 'nx' any of 0, 1, ... l. Once 'nx' has been fixed, we may pick
      // up 'ny' any of: 0, 1, ... (l-nx). Once 'nx' and 'ny' have been fixed, 'nz' can take only
      // the value nz=l-nx-ny. Hence:
      unsigned int pos = 0; // Auxiliar absolute position in the buffer.
      for( unsigned int x=0; x<=order; ++x ){
         for( unsigned int y=0; y<=(unsigned int)((int)order-(int)x); ++y ){
            nx[pos] = x;
            ny[pos] = y;
            nz[pos] = (unsigned int)( (int)order - (int)x - (int)y );
            pos++;
         }
      }
   }
   
   /** This function is intended to compute the multiplicities of each component of a HOT. For a tensor
    of order L, the tensor sumation is defined as terms of the form:
    
    p)          q)         r)             p)            q)            r)
     T_{1,1, ..., 1,2,2, ...,2,3,3, ..., 3} · x·x· ... ·x · y·y· ... ·y · z·z· ... ·z = T_{i_1,i_2,...,i_L} · x^p · y^q · z^r
    
    where p+q+r=L. Since we deal only with hypersymmetric tensors, the term T_{i_1,i_2,...,i_L} is 
    invariant under any permutation of the set of indices {i_1,i_2,...,i_L}:
    
    T_{1,1,2,2,2,3,3} = T_{3,1,2,1,2,3,2} = T_{3,3,2,2,2,1,1} = ...
    
    It is easy to prove that there are exactly (L+1)(L+2)/2 independent terms for a symmetric tensor
    of order L, each of them appearing with a certain multiplicity in the overall sumation. The
    purpose of this method is precisely to compute such factors.
    
    The memory to store the indices has to be allocated externally. Note that HOT form a basis for the
    same functional space as even-order SH, so the memory required can be allocated via:
    
    unsigned int* mu = new unsigned int[ getNumberOfEvenAssociatedLegendrePolynomials(order) ];
    
    The arguments nx, ny, and nz can be computed via the computeSetOfHOTPowers() function.
    */
   void computeMultiplicityOfHOTComponents( const unsigned int order, unsigned long* mu, const unsigned int* nx, const unsigned int* ny, const unsigned int* nz )
   {
      // Get the total number of independent terms of the tensor:
      unsigned int l = (order+1)*(order+2)/2;
      // The multiplicity of a term with powers nx, ny, and nz is easily
      // proven to be: l!/nx!/ny!/nz!. Instead of explicitly computing the
      // factorials, we use the function cumprod:
      for( unsigned int pos=0; pos<l; ++pos ){
         // Find the greater degree among nx, ny, and nz:
         if( nx[pos]>=ny[pos] && nx[pos]>=nz[pos] ){
            mu[pos] = cumprod( nx[pos], order );
            for( unsigned int d=ny[pos]; d>1; --d )
               mu[pos] /= d;
            for( unsigned int d=nz[pos]; d>1; --d )
               mu[pos] /= d;
         }
         else if( ny[pos]>=nx[pos] && ny[pos]>=nz[pos] ){
            mu[pos] = cumprod( ny[pos], order );
            for( unsigned int d=nx[pos]; d>1; --d )
               mu[pos] /= d;
            for( unsigned int d=nz[pos]; d>1; --d )
               mu[pos] /= d;
         }
         else{
            mu[pos] = cumprod( nz[pos], order );
            for( unsigned int d=nx[pos]; d>1; --d )
               mu[pos] /= d;
            for( unsigned int d=ny[pos]; d>1; --d )
               mu[pos] /= d;
         }
      }
   }
   
   /**
    * This function is intended to evaluate a escalar function defined over the unit sphere given by
    * a hypersymmetric higher order tensor. An L-order tensor is evaluated (for 3-D vectors) as:
    *
    *         3      3           3
    *   sum                             A                    x    · x   · ... · x
    *       i_1=1, i_2=1, ..., i_L=1     i_1, i_2, ..., i_L    i_1   i_2         i_L
    *
    * However, since we deal only with hypersymmetric tensors, this summation of 3^L elements can be
    * reduced to a small subset of (L+1)(L+2)/2 independent terms with given muliplicities:
    *
    *         L                 p_j    q_j    r_j
    *    sum      T_j · mu_j · x    · y    · z
    *        j=1
    *
    * The multiplicities mu_j and the powers of x, y, and z for each independent term are obtained
    * via the computeMultiplicityOfHOTComponents() and the computeSetOfHOTPowers() functions,
    * respectively.
    *
    * The meaning of the parameters are as follows:
    *
    *    N:              The number of points where the HOT has to be evaluated
    *    x, y, and z:    Each of them is a vector of length N. For a given position j, 0<=j<N,
    *                    the HOT is evaluated for [x[j],y[j],z[j]]^T, which necessarily has
    *                    norm 1.
    *    val:            The vector with the evaluations of the HOT (output).
    *    order:          The order of the HOT (only even degrees 0, 2, 4, 6, 8, ... are considered)
    *    nx, ny, and nz: The degrees of x, y, and z, respectively, for each free component of the HOT.
    *                    Compute them with the computeSetOfHOTPowers() function.
    *    mu:             The multiplicity of each free component of the tensor. Compute it via the
    *                    computeMultiplicityOfHOTComponents() function.
    *    T:              The actual independent components of the HOT, T_j. Note that we use a
    *                    vnl_matrix to represent it (in fact, an (order+1)·(order+2)/2 x 1 matrix).
    *                    The reason is that these components will be typically obtained from the
    *                    coefficients of an SH expansion with a pass matrix 'M':
    *
    *                            T = M * SH
    *
    *                    where M can be computed with the generateSH2HOTMatrix() function.
    */
   void evaluateHOTFromCoefficients( const unsigned int N, const double* x, const double* y, const double* z, 
                             double* val,
                             const unsigned int order, const unsigned int* nx, const unsigned int* ny, const unsigned int* nz,
                             const unsigned long* mu,
                             const SHMatrixType T )
   {
      // Compute the number of terms of the summation depending on the order
      // of the tensor:
      unsigned int L = (order+1)*(order+2)/2;
      for( unsigned int n=0; n<N; ++n ){ // At each point of the unit sphere
         // Initiallize the output:
         val[n] = 0.000000000000000f;
         // Implement the summation:
         for( unsigned int j=0; j<L; ++j ){ // For each term
            // Initiallize the independent term, x^p·y^q·z^r:
            double cum = 1.000000000000000f;
            // Compute the whole independent term:
            for( unsigned int p=0; p<nx[j]; ++p ) // Implements x^p
               cum *= x[n];
            for( unsigned int q=0; q<ny[j]; ++q ) // Implements y^q
               cum *= y[n];
            for( unsigned int r=0; r<nz[j]; ++r ) // Implements z^r
               cum *= z[n];
            // Sum the result with its corresponding multiplicity:
            val[n] += T[j][0] * mu[j] * cum;
         }
      }
   }
   /** Each term in the sumation of a hypersymmetric, order L tensor product is denoted as:
             p)          q)         r)             p)            q)            r)
    T_{1,1, ..., 1,2,2, ...,2,3,3, ..., 3} · x·x· ... ·x · y·y· ... ·y · z·z· ... ·z = T_{i_1,i_2,...,i_L} · x^p · y^q · z^r,
    
    so that:
   
    T·x^m =  sum_{j=1}^{(L+1)(L+2)/2} T_{i_1(j),i_2(j),...,i_L(j)} · mu(j) · x^p(j) · y^q(j) · z^r(j),
    
    
    where p+q+r=L and mu is the multiplicity of the corresponding term. Since we deal only with 
    hypersymmetric tensors, the term T_{i_1,i_2,...,i_L} is invariant under any permutation of the set
    of indices {i_1,i_2,...,i_L}:
    
    This is the compact form of the sum, which reduces to a far smaller number of terms when compared
    to the standard tensor product:
    
    T·x^m = sum_{i_1=1}^3 sum_{i_2=1}^3 ... sum_{i_L=1}^3 T_{i_1,i_2,...,i_L} · x_{i_1} · x_{i_2} · ... · x_{i_L}
    
    The aim in this function is to find a "reduced tensor" of order L-1 from the original tensor
    of order L. Explicitly, we aim to find the 3x1 vector T·x^{m-1} defined as:
    
     [sum_{i_2=1}^3 sum_{i_2=1}^3 ... sum_{i_L=1}^3 T_{1,i_2,...,i_L} · x_{i_2} · x_{i_2} · ... · x_{i_L}]
     [sum_{i_2=1}^3 sum_{i_2=1}^3 ... sum_{i_L=1}^3 T_{2,i_2,...,i_L} · x_{i_2} · x_{i_2} · ... · x_{i_L}]
     [sum_{i_2=1}^3 sum_{i_2=1}^3 ... sum_{i_L=1}^3 T_{3,i_2,...,i_L} · x_{i_2} · x_{i_2} · ... · x_{i_L}]
    
    This vector is relevant since it represents (or at least is proportional to) the gradient of the 
    tensor product T·x^m with respect to the three free components [x_1,x_2,x_3]^T. At the same time,
    the computation of such gradient is relevant when we aim to find the local extrema of the tensor
    product T·x^m restricted to the unit sphere. Details on this theory can be found in:
       
      Luke Bloy and Ragini Verma,
      "On Computing the Underlying Fiber Directions from the Diffusion Orientation Distribution Function"
      Lecture Notes in Computer Science 5241: Medical Image Computing and Computer Assisted Intervention
      (Part I) pp. 1-8. New York, September 2008. Eds: D. Metaxas et al.
    
    Note that, since we obtain tensors of order L-1, this function IS NOT restricted to even order tensors.
    It returns a vnl_matrix of size L(L+1)/2 x 1 (the free components of the reduced tensor of order L-1)
    whose elements correspond to the free components of the reduced tensor extracted from the original
    tensor; to compute the powers of x_i and the corresponding multiplicities for each of them, the 
    functions computeSetOfHOTPowers() and computeMultiplicityOfHOTComponents() should be used. This
    function has to be subsequently called for each of the three components of T·x^{m-1}. For example:
    
      SHMatrixType Tx = computeReducedTensor( T, order, 0 );
      SHMatrixType Ty = computeReducedTensor( T, order, 1 );
      SHMatrixType Tz = computeReducedTensor( T, order, 2 );
    
    For all three reduced tensors, the powers and multiplicities are obtained as:
    
      computeSetOfHOTPowers( order-1, nx, ny, nz );
      computeMultiplicityOfHOTComponents( order-1, mu, nx, ny, nz );
    
    where the vectors have lengths:
    
      unsigned int* nx = new unsigned int[order*(order+1)/2];
    
    */
   SHMatrixType computeReducedTensor( const SHMatrixType original, const unsigned int order, const unsigned int coordinate )
   {
      // In practice, since the tensors we are dealing with are all shyprsymmetric, the algorithm may
      // be kept quite simple: for the component selected, we keep the first index of the tensor equal
      // to the coordinate chosen and vary all others. This implies that we only need to remove the
      // L+1 free components of the original tensor for whic the power of the desired coordinate is 0.
      // To prove that this is consistent, simply check that the original tensor has (L+1)(L+2)/2 free
      // components and the reduced tensor has L(L+1)/2, so we loose (L+1)(L+2)/2 - L(L+1)/2 = L+1
      // free components (those which imply that the power of the desired component is 0. In fact, the
      // problem is as simple as copying the elements of the old tensor in the new one in the same
      // order, and remove those which should not be included in the reduced HOT.
      SHMatrixType reduced( order*(order+1)/2, 1 ); // Initiallize the new tensor
      unsigned int pos_n = 0; // Auxiliar absolute position in the buffer of the new tensor.
      unsigned int pos_o = 0; // Auxiliar absolute position in the buffer of the old tensor
      unsigned int powers[3]; // Auxiliar vector to check if the component has to be included 
      for( unsigned int x=0; x<=order; ++x ){
         for( unsigned int y=0; y<=(unsigned int)((int)order-(int)x); ++y ){
            powers[0] = x;
            powers[1] = y;
            powers[2] = (unsigned int)( (int)order - (int)x - (int)y );
            // The value is included in case the power associated to the component
            // chosen is not 0.
            if( powers[coordinate] != 0 ) 
               reduced(pos_n++,0) = original(pos_o,0);
            ++pos_o;
         }
      }
      return reduced;
   }
   
   
   // /* THESE METHODS DO NOT COMPILE WITH SLICER BECAUSE OF VNL_RPOLY_ROOTS SYMBOL NOT FOUND. WE DO NOT NEED THEM RIGHT NOW...
   
   /** The purpose of this function is to find the local extrema of a function defined over the unit sphere by
    means of a Higher Order Tensor. When using Spherical Harmonics, the problem "reduces" to differentiate the
    SH expansion and then solve for theta and phi. However, this problem is not computationally tractable.
    Instead, the solution introduce in:
    
      Luke Bloy and Ragini Verma,
      "On Computing the Underlying Fiber Directions from the Diffusion Orientation Distribution Function"
      Lecture Notes in Computer Science 5241: Medical Image Computing and Computer Assisted Intervention
      (Part I) pp. 1-8. New York, September 2008. Eds: D. Metaxas et al.
    
    is to convert the SH expansion to an equivalent HOT and use the method of Lagrange multipliers, giving
    rise to the vector equation (with three components):
    
      T·e^{l-1} = lambda·e, s.t. e^t·e = 1
    
    (see the help on the function "computeReducedTensor" for further details on the notation). Vectors 
    e=[x,y,z]^t and scalars lambda fulfilling the previous equation are respectively called Z-eigenvectors
    and Z_eignevalues in:
    
      L. Qi, "Eigenvalues of a real supersymmetric tensor"
      Journal of Symbolic Computation 40(6): 1302-1324. 2005
    
    Details on the theory behind this function can be found in the two references above. Like in
    Bloy & Verma (2008), we look for solutions of one of three forms: (1,0,0), (t,1,0), and (t,u,1),
    being able to reduce a system of four equations with four unknowns to (at most) a system of three
    equations with three unknowns.
    */
   /*
   void computeLocalExtrema( const unsigned int order, const SHMatrixType hot, std::vector< vnl_vector<double> >& output )
   {
      // First of all, from the HOT, we have to find the three reduced tensors corresponding to each
      // of the three components of T·x^{l-1}. Each of them has order order-1:
      SHMatrixType rtx = computeReducedTensor( hot, order, 0 );
      SHMatrixType rty = computeReducedTensor( hot, order, 1 );
      SHMatrixType rtz = computeReducedTensor( hot, order, 2 );
      // Note that eq. (4) in Bloy et al. (2008) is given for the extended form of the HOT, but we
      // represent HOT only with their independent terms and multiplicities. Hence, we need to
      // compute them in order to arrange the polynomial equations:
      unsigned int*  nx = new unsigned int[order*(order+1)/2];
      unsigned int*  ny = new unsigned int[order*(order+1)/2];
      unsigned int*  nz = new unsigned int[order*(order+1)/2];
      unsigned long* mu = new unsigned long[order*(order+1)/2];
      computeSetOfHOTPowers( order-1, nx, ny, nz );
      computeMultiplicityOfHOTComponents( order-1, mu, nx, ny, nz );
      // The degrees of the polynomials are no longer necessary:
      delete[] nx;
      delete[] ny;
      delete[] nz;
      // Sequentially call the solvers for each kind of solutions:
      output.clear();
      computeSolutionForm100( order, rtx, rty, rtz, output );
      computeSolutionFormt10( order, rtx, rty, rtz, output, mu );
      computeSolutionFormtu1( order, rtx, rty, rtz, output, mu );
      // The multiplicities are no longer necessary:
      delete[] mu;
      return;
   }
   */
   /** See th help on "computeLocalExtrema". This function is used to find Z-eigenvectors of the form [1,0,0]^T
       IMPORTANT: "order" refers to the order of the original tensor, and hence hx, hy, hz are HOT of order "order-1"
    */
   /*
   void computeSolutionForm100( const unsigned int order, const SHMatrixType hotx, const SHMatrixType hoty, const SHMatrixType hotz,
                         std::vector< vnl_vector<double> >& output )
    
   {
      // The only unknown in this case is the Z-eigenvalue, lambda. Note that y=z=0, so most of the terms
      // in the summation of the tensor are simply 0. In fact, only the last term in the hotx subtensor, 
      // the one corresponing to x^{order-1}, is not null. Hence, the equation reduces to:
      //
      //     T_{0,0,...} = lambda·1,
      //
      // which has the trivial solution lambda = T_{0,0,...}. However, for this solution to be valid we need
      // that the remaining two equations are fulfilled, i.e.:
      //
      //     [T·e^{order-1}]_y = 0, [T·e^{order-1}]_z = 0
      //
      // Due to the very simple form of 'e' in this case, it is very simple to check this condition: it is
      // enough to test if the last components of hoty and hotz are both zero.
      if( fabs(hoty(order*(order+1)/2-1,0))<1e-12 && fabs(hotz(order*(order+1)/2-1,0))<1e-12 ){
         vnl_vector<double> sol(3);
         sol[0] = 1.000000000000000f;
         sol[1] = 0.000000000000000f;
         sol[2] = 0.000000000000000f;
         output.push_back( sol );
      }
      return;
   }
   */
   /** See th help on "computeLocalExtrema". This function is used to find Z-eigenvectors of the form [t,1,0]^T
    IMPORTANT: "order" refers to the order of the original tensor, and hence hx, hy, hz are HOT of order "order-1"
    */
   /*
   void computeSolutionFormt10( const unsigned int order, const SHMatrixType hotx, const SHMatrixType hoty, const SHMatrixType hotz,
                         std::vector< vnl_vector<double> >& output, const unsigned long* mu )
   {
      // Now, we have two unknowns: the Z-eignevalue 'lambda' and the parameter 't'. Hence, we need two equations
      // (those corresponding to [T·e^{order-1}]_x and [T·e^{order-1}]_y and additionally we need to check if the
      // third equation on [T·e^{order-1}]_z is fulfilled. Since we assume that z=0, only those terms in the summations 
      // comprising powers only of x and y will be present. For order "order-1", it is trivial to show that there are
      // exactly "order" such terms. On the other hand, the unknown lambda can be easily eliminated:
      //
      //    [T·e^{order-1}]_x = lambda·x
      //    [T·e^{order-1}]_y = lambda·y
      //
      //       then
      //
      //    [T·e^{order-1}]_x / x = [T·e^{order-1}]_y / y, or
      //    [T·e^{order-1}]_x·y - [T·e^{order-1}]_y·x = 0 => [T·e^{order-1}]_x - [T·e^{order-1}]_y·x = 0
      //
      // and we have a unique polynomial in the unknown 'x' ('y' is assumed to be 0). Note that this polynomial
      // has 2·order terms. Hence:
      //
      // The maximum degree of the polynomial will be given by [T·e^{order-1}]_y·x; since the maximum degree of
      // [T·e^{order-1}]_y in 'x' is order-1, the global degree is order. Hence, we need order+1 coefficients:
      vnl_vector<double> coefficients(order+1);
      coefficients.fill(0.000000000000000f);
      // The coefficients are stored in decreasing powers; we fill the vector:
      unsigned int pos = 0; // Auxiliar absolute position in the HOT.
      for( unsigned int x=0; x<=order-1; ++x ){ // For all possible values of x
         for( unsigned int y=0; y<=(unsigned int)((int)order-1-(int)x); ++y, ++pos ){
            unsigned int z = (unsigned int)( (int)order - 1 - (int)x - (int)y );
            // We have to include this term only in case the power of the z
            // component is 0, so that the corresponding term is not null:
            if( z==0 ){
               coefficients[order-x]   += hotx(pos,0); // corresponding term in [T·e^{order-1}]_x
               coefficients[order-x-1] -= hoty(pos,0); // corresponding term in -[T·e^{order-1}]_y·x
            }
         }
      }
      // Make sure the leading coefficient is not 0 (not allowed by vnl_rpoly_roots):
      for( pos=0; pos<=order; ++pos ){
         if( fabs(coefficients[pos])>1e-12 )
            break;
      }
      if( pos>order )
         return;
      else if( pos>0 )
         coefficients = coefficients.extract( (unsigned int)((int)order+1-(int)pos), pos );
      // Now, we can call the vnl implementation to find the roots of the polynomial in x:
      vnl_rpoly_roots solver( coefficients ); // The constructor already does all computations.
      // We are only interested in real solutions:
      vnl_vector<double> reals = solver.realroots();
      // Extract and normalize the solutions:
      vnl_vector<double> sol(3);
      for( unsigned int k=0; k<reals.size(); ++k ){
         sol[0] = reals[k];
         sol[1] = 1.000000000000000f; // Remember that 'y' and 'z' are both fixed; y=1
         sol[2] = 0.000000000000000f; // Remember that 'y' and 'z' are both fixed; z=0
         // For this solution to be valid, remember that it must fulfill the last
         // equation, the one in [T·e^{order-1}]_z. We may test the "order" non-null
         // terms in hotz to check if the sum is in fact zero:
         pos          = 0;
         double sum   = 0.000000000000000f; // Auxiliar value to compute the tensor sum
         double power = 1.000000000000000f; // Auxiliar value to compute x^p(r) (see below)
         for( unsigned int x=0; x<=order-1; ++x ){ // For all possible values of x
            for( unsigned int y=0; y<=(unsigned int)((int)order-1-(int)x); ++y, ++pos ){
               unsigned int z = (unsigned int)( (int)order - 1 - (int)x - (int)y );
               if( z==0 ){
                  // Compute T_r·mu(r)·x^p(r)·y^q(r) (note that the power of 'z' is 0)
                  // Since 'y' is assumed to be 1, it reduces to T_r·mu(r)·x^p(r); note
                  // that x^p(r) is computed externally
                  sum += hotz(pos,0)*mu[pos]*power;
               }
            }
            power *= sol[0];
         }
         // If the last equation is fulfilled, we have a correct solution:
         if( fabs(sum)<1e-12 ){
            // Normalize to achieve norm 1:
            sol /= (sol.magnitude());
            output.push_back( sol );
         }
      }
      return;
   }
    */
   
   /** See th help on "computeLocalExtrema". This function is used to find Z-eigenvectors of the form [t,u,1]^T
    IMPORTANT: "order" refers to the order of the original tensor, and hence hx, hy, hz are HOT of order "order-1"
    */
   /*
   void computeSolutionFormtu1( const unsigned int order, const SHMatrixType hotx, const SHMatrixType hoty, const SHMatrixType hotz,
                         std::vector< vnl_vector<double> >& output, const unsigned long* mu )
   {
      // In this case we have three unknowns: 'x', 'y', and the Z-eignevalue 'lambda', while 'z' is assumed to be 1.
      // This is the most general case, and will be also the most frequent. We have to consider all the terms in the
      // HOT summation, since none of them will be null a priori. As in the previous case, we can eliminate lambda:
      //    [T·e^{order-1}]_x = lambda·x
      //    [T·e^{order-1}]_y = lambda·y
      //    [T·e^{order-1}]_z = lambda·z
      //
      //       then
      //
      //    [T·e^{order-1}]_x·z = [T·e^{order-1}]_z·x,
      //    [T·e^{order-1}]_y·z = [T·e^{order-1}]_z·y. 
      //
      //       finally
      //
      //    [T·e^{order-1}]_x - [T·e^{order-1}]_z·x = 0, (1)
      //    [T·e^{order-1}]_y - [T·e^{order-1}]_z·y = 0. (2)
      //
      // For each equation of the system, we have the order·(order+1)/2 terms of the left hand side and
      // the corresponding order·(order+1)/2 terms from the right hand side, order·(order+1) as a total:
      vnl_matrix<unsigned int> p1( order*(order+1), 2 );  // The degrees for equation (1)
      vnl_matrix<unsigned int> p2( order*(order+1), 2 );  // The degrees for equation (2)
      vnl_vector<double>       c1( order*(order+1) );     // The coefficients for equation (1)
      vnl_vector<double>       c2( order*(order+1) );     // The coefficients for equation (2)
      // Now, proceed to fill the coefficients and the powers:
      unsigned int pos    = 0;                 // Auxiliar absolute position in the HOT.
      unsigned int offset = order*(order+1)/2; // Offset to access the right hand side terms
      for( unsigned int x=0; x<=order-1; ++x ){ // For all possible values of x
         for( unsigned int y=0; y<=(unsigned int)((int)order-1-(int)x); ++y, ++pos ){
            p1(pos,0) = p2(pos,0) = x;   // first unknown has degree 'x'
            p1(pos,1) = p2(pos,1) = y;   // second unknown has degree 'y'
            p1(pos+offset,0)      = x+1; // the right hand side of p1 is smultiplied by 'x', hence we add 1
            p1(pos+offset,1)      = y;   // the degree in the right hand side term is the same
            p2(pos+offset,0)      = x;   // the degree in the right hand side term is the same
            p2(pos+offset,1)      = y+1; // the right hand side of p2 is smultiplied by 'y', hence we add 1
            c1[pos]               =  hotx(pos,0) * mu[pos];
            c2[pos]               =  hoty(pos,0) * mu[pos];
            c1[pos+offset]        = -hotz(pos,0) * mu[pos];
            c2[pos+offset]        = -hotz(pos,0) * mu[pos];
         }
      }
      // Now, we can arrange the polynomial equations and the system of equations:
      vnl_real_npolynomial eq1( c1, p1 );
      vnl_real_npolynomial eq2( c2, p2 );
      vcl_vector<vnl_real_npolynomial*> eqsystem( 2 );
      eqsystem[0] = &eq1;
      eqsystem[1] = &eq2;
      // Now, we can call the vnl implementation to solve the system for the unknowns x, y, and z:
      vnl_rnpoly_solve solver( eqsystem ); // The constructor already does all computations.
      // We are only interested in real solutions;
      vcl_vector<vnl_vector<double>*> reals = solver.realroots();
      // Extract and normalize the solutions:
      vnl_vector<double> sol(3);
      vcl_vector<vnl_vector<double>*>::iterator rit;
      for( rit=reals.begin(); rit!=reals.end(); ++rit ){
         sol[0] = (*rit)->get(0);
         sol[1] = (*rit)->get(1);
         sol[2] = 1.000000000000000f; // Remember that 'z' is 1
         // This solution is always valid, so no additional checking is performed:
         sol /= (sol.magnitude());
         output.push_back( sol );
      }
      return;
   }
    */
   
} // End namespace shmaths

#endif // #ifndef _sh2hot_cxx
