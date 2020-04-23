1==1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% THE STARTING POINT FOR THIS FILE WAS HANNAH's functionsHex.m found in    %%
%%   FILES\private_notes\2016-07-04-hannah_thesis_codes                     %%
%%       \Spin_Clusters\Spin-Clusters-master                                %%
%%       \ExampleHexagon\Fig6_16\functionsHex.m                             %%
%% For an eplanation of those codes see                                     %%
%%   FILES\private_notes\2016-07-04-hannah_thesis_codes                     %%
%%       \Spin_Clusters\Spin-Clusters-master                                %%
%%       \ExampleHexagon\ReadMe.txt                                         %%
%% In this file I have added the follwoing functions:
%%
%%   gen_pic_mat                                                            %%
%%   gen_r_vec                                                              %%
%% JQ, 2 September 2019                                                     %%
%% IMPORTANT NOTE: added 3 September 2019
%%
%% As it stands, gen_r_vec(N) gives the position vectors used by Hannah
%% only for N=6. For N=10, it gives the same position vectors but with
%% different indices. For N=4,8 the molecule appears rotated. This can be
%% seen by comparing r_vec_comparison_HANNAH_v01.pdf to
%% . Because of this, the code currently gives
%% correct results ONLY FOR N=4.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wrapper function for gen_pic_mat wbhich takes the parameters in a single
% vector (for use with function gen_A in functions_PCA_v01-JQ.m):
%
function pic_mat=gen_pic_mat_wrapped(params)
    N=params(1);
    Gamma=params(2);
    Delta=params(3);
    h=params(4);
    T=params(5);
    q_max=params(6);
    Nk=params(7);
    pic_mat=gen_pic_mat(N,Gamma,Delta,h,T,q_max,Nk);
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gen_pic_mat
%
%   Calculation of the total diffuse magnetic scattering function S(q)
%   off a planar molecule with N spins having nearest-neighbour anisotropic
%   Heisenberg magnetic interactions at temperature T and in a magnetic
%   field h perpendicular to the plan of the molecule.
%
%   The scattering vector q is assumed to be in the plane of the molecule.
%
% INPUT:
%
%   r_vec   N x 2 matrix giving the positions of all N atoms
%           forming the molecule (the atoms are assumed to be on
%           the xy plane and their positions are given in Cartesian
%           coordinates)
%
%   Gamma   Real number giving the in-plane anisotropy of the spin-spin
%           interactions (for Delta = 0, Gamma=0 is the XY model,
%           Gamma=1 is the Ising model)
%
%   Delta   Real number giving the out-of-plane anisotropy of the
%           spin-spin interactions (for Gamma = 0, Delta = 0 is the XY model,
%           Delta = 1 is the Ising model)
%
%   h       Applied magnetic field, as a real number (field assumed to be
%           in the z direction)
%
%   T       Temperature, as a real number
%
%   Nk      Parameters giving the values taken by the scattering vector
%   q_max   q=[qx,qy]. Each component qx,qy takes 2*Nk uniformly-spaced
%           values between -q_max and q_max, inclusive (this automaitcally
%           avoids q=[0,0], whose corresponding unit vector is not
%           well-defined).
%
% OUTPUT:
%
%   pic_mat Nk x Nk matrix giving the values of the scattering function S(q)
%           at each value of q. Note the values of q are not stored in the
%           matrix - only the values of S(q) are.
%
% USAGE EXAMPLE:
%
%   octave:11> source "functionsHex_v02-JQ.m";
%   octave:12> N=6;Gamma=0.5;Delta=0.0;h=0.0;T=0.1;q_max=2.5*pi;Nk=30;
%   octave:13> pic_mat=gen_pic_mat(gen_r_vec(N),Gamma,Delta,h,T,q_max,Nk);
%   octave:14> save "test_03.dat" pic_mat
%
%   gnuplot> set pm3d; unset surface; set view map; set size ratio 1; unset key;
%   gnuplot> Nk=30;
%   gnuplot> set xtics ("-q_{max}a" 0,"0" (2*Nk-1)/2,"q_{max}a" 2*Nk-1);
%   gnuplot> set ytics ("-q_{max}a" 0,"0" (2*Nk-1)/2,"q_{max}a" 2*Nk-1);
%   gnuplot> set title "(a = distance between n.n. atoms within the molecule)";
%   gnuplot> splot [0:2*Nk-1][0:2*Nk-1] "test_03.dat" matrix
%
function pic_mat=gen_pic_mat(N,Gamma,Delta,h,T,q_max,Nk)
%
    r_vec=gen_r_vec(N);
    N = length(r_vec);
    dk=2*q_max/(2*Nk-1);
    qx=(-q_max:dk:q_max);
    qy=qx;
    sigmaz=[1,0;0,-1];
%
   [states_array,energies_vec]=DiagpbHmat(N,Gamma,Delta,h);
   partZ=Z(T,energies_vec);
%
   SQxx=zeros(length(qx),length(qy));
   SQyy=zeros(length(qx),length(qy));
   SQzz=zeros(length(qx),length(qy));
   SQxy=zeros(length(qx),length(qy));
   SQyx=zeros(length(qx),length(qy));
%
       for delta=1:N
          for deltap=1:N
              thetad=theta(delta,N);
              thetadp=theta(deltap,N);
              sigmax1=sigmaxp(thetad);
              sigmay1=sigmayp(thetad);
              sigmax2=sigmaxp(thetadp);
              sigmay2=sigmayp(thetadp);
              FCorrxx(delta,deltap)=finiteCorr(delta,deltap,sigmax1,sigmax2,N,states_array,energies_vec,T,partZ);
              FCorryy(delta,deltap)=finiteCorr(delta,deltap,sigmay1,sigmay2,N,states_array,energies_vec,T,partZ);
              FCorrzz(delta,deltap)=finiteCorr(delta,deltap,sigmaz,sigmaz,N,states_array,energies_vec,T,partZ);
              FCorrxy(delta,deltap)=finiteCorr(delta,deltap,sigmax1,sigmay2,N,states_array,energies_vec,T,partZ);
              FCorryx(delta,deltap)=finiteCorr(delta,deltap,sigmay1,sigmax2,N,states_array,energies_vec,T,partZ);
         endfor
      endfor
%
   for xx=1:length(qx)
       for yy=1:length(qy)
           Sqxx=0;
           Sqyy=0;
           Sqzz=0;
           Sqxy=0;
           Sqyx=0;
           for delta=1:N
              for deltap=1:N
                  vec=[qx(xx),qy(yy)]*(r_vec(delta,:)-r_vec(deltap,:))';
                  Sqxx=Sqxx+(FCorrxx(delta,deltap))*exp(i*vec);
                  Sqyy=Sqyy+(FCorryy(delta,deltap))*exp(i*vec);
                  Sqzz=Sqzz+(FCorrzz(delta,deltap))*exp(i*vec);
                  Sqxy=Sqxy+(FCorrxy(delta,deltap))*exp(i*vec);
                  Sqyx=Sqyx+(FCorryx(delta,deltap))*exp(i*vec);
                  qxhat=xhat(qx(xx),qy(yy));
                  qyhat=xhat(qy(yy),qx(xx));
                  SQxx(yy,xx)=(1-qxhat^2)*Sqxx;
                  SQyy(yy,xx)=(1-qyhat^2)*Sqyy;
                  SQzz(yy,xx)=Sqzz;
                  SQxy(yy,xx)=-(qxhat*qyhat)*Sqxy;
                  SQyx(yy,xx)=-(qyhat*qxhat)*Sqyx;
              endfor
           endfor
       endfor
   endfor
%
   pic_mat=real(SQxx)+real(SQyy)+real(SQzz)+real(SQxy)+real(SQyx);
%
endfunction

% Position vectors of atoms in molecule with N atoms,
% expressed in units of the distance a between n.n. atoms.
%
function r_vec=gen_r_vec(N)
    % Distance between n.n. atoms in units of the molecule's radius:
    a=2*sin(pi/N);
    % Angle increment from one atom to the next:
    dphi=2*pi/N;
    % Angles at which all N atoms are found:
    phi=(pi/2-0.5*dphi:-dphi:-3*pi/2+0.5*dphi)';
    % Coordinates of all atoms (in units of the radius):
    x=cos(phi);
    y=sin(phi);
    r_vec=[x,y];
    % Divide by a to obtain final result in the desired units:
    r_vec=r_vec/a;
endfunction

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%
% FUNCTION: PBspinConstruction : same as above for PERIODIC BOUNDARY CONDITIONS
%           This function creates spin matrices that describe the interaction and field
%           componants in the Hamiltonian
% INPUT:... sigma     Pauli Matrix (2x2). for x y z
%           N         Positive Interger. Size of chain. has limiting factors. N=13 is max.
%
% OUTPUT:...Spin Matrix (2^N,2^N) for S(alpa)i*S(alpha)i+1
%           Field Matrix (2^n,2^n) transverse field componant.
%%%%%%%%
%%%%%%%%%
function [pbSpin_matrix,F_matrix]=PBspinConstruct(sigma,N)

Spin_matrix=zeros((2^N),(2^N));

for ii=1:N-1

	interactionMat=kron(eye(2^(ii-1),2^(ii-1)),(kron(sigma,eye(2^(N-ii),2^(N-ii)))))*kron(eye(2^(ii),2^(ii)),(kron(sigma,eye(2^(N-ii-1),2^(N-ii-1)))));

	Spin_matrix=Spin_matrix+interactionMat;

endfor
  pb=kron(eye(2^(N-1),2^(N-1)),sigma)*kron(sigma,eye(2^(N-1),2^(N-1)));
  pbSpin_matrix=Spin_matrix+pb;

% field componant. NOTE: this assumes a transverse field in the Z direction

   F_matrix=zeros(2^N,2^N);

   for jj=1:N

	    fieldElement=kron(eye(2^(jj-1)),(kron(([1,0;0,-1]),eye(2^(N-jj),2^(N-jj)))));
	    F_matrix=F_matrix+fieldElement;

   endfor

end
%%%%%%%%%
%%%%%%%%%%
% FUNCTION: pbHmat: same as above for PERIODIC BOUNDARY CONDITIONS
%          This function calls in function spinConstruct threee times for SxSx, SySy and SzSz.
%          This means that it can be applied to the xyx model also. It takes the spin
%          interaction matrices and field componant and evaluates them in the Heisenberg
%          Hamiltonian.
% INPUT:...N      Positive Integer. Size of finite Chain.
%          Gamma  Positive number 0<Gamma<1. Anisotropy for x and y
%          Delta  Positive number. Anisotropy for z. taken as zero for xy model.
%          h      Positive number. Effective field. i.e h=1 is the critical field.
% OUTPUT...The Hamiltonian matrix (2^N,2^N)
function Hmatrix=pbHmat(N,Gamma,Delta,h)
% combime for the Hamiltonian. INPUT: Gamma, Delta and field call in above function, require input N

   [SxSx,F_matrix]=PBspinConstruct([0,1;1,0],N);   % call in for x,y,z
   [SySy,F_matrix]=PBspinConstruct([0,i;-i,0],N);
   [SzSz,F_matrix]=PBspinConstruct([1,0;0,-1],N);

   Vx=SxSx*(1+Gamma);     % add anisotropy for x,y and z
   Vy=SySy*(1-Gamma);
   Vz=SzSz*Delta;

   Hmatrix=(1/4)*(Vx+Vy+Vz)-(1/2)*h*F_matrix;       % construct hamiltonian matrix with field.

end
%%%%%%%%%%
%%%%%%%%%%
% FUNCTION: DiagpbHmat: same as above for PERIODIC BOUNDARY CONDITIONS
%           This calls in Hmat so would be the only required function to call by a script to
%           generate eigenvalues and eigen vecors for a system for a given range of h. h can be
%           a singular number or a vector. It does this through exact diagonilisation
% INPUT:... N      Positive Integer. Size of Finite chain
%           Gamma  Positive number. 0<Gamma<1. Anisotropy for x and y
%           Delta  Positive number. Anisotropy for z. taken as zero for xy model.
%           h_vec  Positive number or vector. Effective field. i.e h=1 is the critical field.
% OUTPUT... eVec   Matrix (2^N,2^N) of eigen vectors.
%           eVal   Vector 2^N of eigen values.
function [eVec,eVal]=DiagpbHmat(N,Gamma,Delta,h_vec)
for ii=1:length(h_vec) % h_vec can be a single field value of a vector range of them.

	h=h_vec(ii); %(not 100% on this one for h as a vector) check this
	[e_vector,e_value]=eig(pbHmat(N,Gamma,Delta,h));
	eVec(:,:,ii)=e_vector;
	eVal(:,ii)=diag(e_value);

endfor
end
%%%%%%%%%%
%%%%%%%%%%
% FUNCTION: theta: rotation angle per site for a plaquette
% INPUT:    delta   site delta 1:4
% OUTPUT:   thetan  theta for delta
%%%%%%%%%%
%%%%%%%%%%
function thetan=theta(delta,N)
   thetan=((2*delta)-1)*(pi/N);
endfunction

%%%%%%%%%%
%%%%%%%%%%
% FUNCTION: sigmax: calculates x prime sigma x per site dependent on the rotation
% INPUT:    theta in itself dependent on site
% OUTPUT:   rotated operator
%%%%%%%%%%
%%%%%%%%%%
function sigmax=sigmaxp(theta)

px=[0,1;1,0];
py=[0,-i;i,0];

Theta=(2*pi)-theta;

sigmax=px*cos(Theta)-py*sin(Theta);

end
%%%%%%%%%%
%%%%%%%%%%
% FUNCTION: sigmay: calculates y prime sigma y per site dependent on the rotation
% INPUT:    theta in itself dependent on site
% OUTPUT:   rotated operator
%%%%%%%%%%
%%%%%%%%%%
function sigmay=sigmayp(theta)

px=[0,1;1,0];
py=[0,-i;i,0];

Theta=(2*pi)-theta;

sigmay=px*sin(Theta)+py*cos(Theta);

end
%%%%%%%%%%
%%%%%%%%%%
% FUNCTION: corrPair
%           Constructs a realspace correlator operator for a given pair in a finite system.
%           To be called in by pureCorrPoint
% INPUT:... delta   Positive Integer. Position i
%           deltap  Positive Integer. Position j
%           sigma_1 and sigma_2 pauli matrices. what type of correlator do you want? xx, xy etc
%           N       Positive Integer. size of Finite system.
% OUTPUT... correlatorPair   matrix (2^N,2^N) for a set correlation between i and j sites.
function correlatorPair=corrPair(delta,deltap,sigma_1,sigma_2,N)

corr_1=zeros((2^N),(2^N));   % initiate
corr_2=zeros((2^N),(2^N));


      corr_1(:,:)=kron(eye(2^(delta-1),2^(delta-1)),(kron(0.5*sigma_1,eye(2^(N-delta),2^(N-delta)))));        % position i

      corr_2(:,:)=kron(eye(2^(deltap-1),2^(deltap-1)),(kron(0.5*sigma_2,eye(2^(N-deltap),2^(N-deltap)))));       % position j

         correlatorPair(:,:)=corr_1(:,:)*corr_2(:,:);

end

%%%%%%%%%%%
%%%%%%%%%%%
% FUNCTION: pureCorr
%           Evaluates the realspace correlation value for a given pair of sites i and j.
%           assumes the system is in a pure state. to then be fed into finiteCorr for finite
%           temperatture results. Groundstate given here for T=0 when using states_array(:1) %           only. Calls in operator function (above).
% INPUT:... delta   Positive Integer. Position i
%           deltap  Positive Integer. Position j
%           sigma_1 and sigma_2 pauli matrices. what type of correlator do you want? xx, xy etc
%           N       Positive Integer. size of Finite system.
%           states_array. matrix (2^N,2^N) are the eigenstates calculated from diagonalising
%           the Hamiltonian. Using function DiagHmat.
%           NOTE: must take states_array as a matrix, can't take one vector at a time. this
%           feature is integrated into finiteCorr function as to not require this one.
% OUTPUT:...pureCorrPoint   vector (2^N) correlator value per state for a given pair of sites.
function pureCorrPoint=pureCorr(delta,deltap,sigma_1,sigma_2,N,states_array)

   correlatorPair=corrPair(delta,deltap,sigma_1,sigma_2,N);   % compute operator pair
                                           % compute expectation value for each state

      for ii=1:2**N
               pureCorrPoint(ii,:)=(conj(states_array(:,ii)'))*(correlatorPair*states_array(:,ii));
      endfor

end
%%%%%%%%%%%
% FUNCTION: Z
%           Evaluates the Partition Function for a given set of eigen values.
% INPUT:... T               Positive number. Temperature as a percentage of J.
%           energies_array  vector 2^N.      Eigenvalues of the system.
% OUTPUT:...partitionFunction Positive number.
function partitionFunction=Z(T,energies_array)

   [x,y] = size (energies_array);

        partitionFunction = 0;
        for ii=1:x
                partitionFunction = partitionFunction + exp(-(energies_array(ii)/T));
        endfor


end
%%%%%%%%%%%%
%%%%%%%%%%%%
% FUNCTION: finiteCorr
%           Evaluates the correlators in pair for a given temperature T. calls in pureCorr.
% INPUT:... delta   Positive Integer. Position i
%           deltap  Positive Integer. Position j
%           sigma_1 and sigma_2 pauli matrices. what type of correlator do you want? xx, xy etc
%           N       Positive Integer. size of Finite system.
%           states_array. matrix (2^N,2^N) are the eigenstates calculated from diagonalising
%           the Hamiltonian. Using function DiagHmat.
%           T               Positive number. Temperature as a percentage of J.
%           energies_array  vector 2^N.      Eigenvalues of the system.
%           partZ   Positive number. Partition Function
% OUTPUT:...finiteCorrelator. A number. the realspace correlation value for a given pair difined
%           by delta and deltap position at a given temperature for teh original set of
%           parameters. Gamma, Delta, field etc.
function finiteCorrelator=finiteCorr(delta,deltap,sigma_1,sigma_2,N,states_array,energies_vec,T,partZ)


 correlatorPair=corrPair(delta,deltap,sigma_1,sigma_2,N);
%[x,y]=size(energies_vec);

finCorr=0;
	for ii=1:2**N
	   pureCorrPoint(ii,:)=(conj(states_array(:,ii)'))*(correlatorPair*states_array(:,ii));
       value=(pureCorrPoint(ii)*exp((-energies_vec(ii))/T))/partZ;
       finCorr=finCorr+value;

	endfor
   finiteCorrelator=finCorr;
endfunction

% FUNCTION: gsCorr
%           Evaluates the correlators in pair for a given temperature T. calls in pureCorr.
% INPUT:... delta   Positive Integer. Position i
%           deltap  Positive Integer. Position j
%           sigma_1 and sigma_2 pauli matrices. what type of correlator do you want? xx, xy etc
%           N       Positive Integer. size of Finite system.
%           state   Vector (2**N) the groundstate
% OUTPUT:...gsCorrelator. Groundstate correlator value for pair of delta deltap.
function gsCorrelator=gsCorr(delta,deltap,sigma_1,sigma_2,N,state)

   correlatorPair=corrPair(delta,deltap,sigma_1,sigma_2,N);

   gsCorrelator=(conj(state'))*(correlatorPair*state);

end


%%%%%%%%%%%
%%%%%%%%%%%
%
% FUNCTION xhat
%          unit vector componant
% INPUT    qx qy
%
% Output   qhat
%%%%%%%%%%%
function qhat=xhat(q1,q2)

qhat=q1/(sqrt((q1^2)+(q2^2)));


endfunction
