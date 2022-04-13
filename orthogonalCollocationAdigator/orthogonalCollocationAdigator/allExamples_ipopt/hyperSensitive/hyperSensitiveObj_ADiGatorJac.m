% This code was generated using ADiGator version 1.4
% ©2010-2014 Matthew J. Weinstein and Anil V. Rao
% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ 
% Contact: mweinstein@ufl.edu
% Bugs/suggestions may be reported to the sourceforge forums
%                    DISCLAIMER
% ADiGator is a general-purpose software distributed under the GNU General
% Public License version 3.0. While the software is distributed with the
% hope that it will be useful, both the software and generated code are
% provided 'AS IS' with NO WARRANTIES OF ANY KIND and no merchantability
% or fitness for any purpose or application.

function obj = hyperSensitiveObj_ADiGatorJac(z)
global ADiGator_hyperSensitiveObj_ADiGatorJac
if isempty(ADiGator_hyperSensitiveObj_ADiGatorJac); ADiGator_LoadData(); end
Gator1Data = ADiGator_hyperSensitiveObj_ADiGatorJac.hyperSensitiveObj_ADiGatorJac.Gator1Data;
% ADiGator Start Derivative Computations
%User Line: % Computes the objective function
global psStuff nstates ncontrols 
%User Line: global
%User Line: % Extract the number of LGR points, the LGR differentiation matrix,
%User Line: % the LGR points, and the LGR weights
D = psStuff.D;
%User Line: D = psStuff.D;
tau = psStuff.tau;
%User Line: tau = psStuff.tau;
w = psStuff.w;
%User Line: w = psStuff.w;
%User Line: % Decompose the NLP decision vector into pieces containing
%User Line: %  - the state
%User Line: %  - the control
%User Line: %  - the initial time
%User Line: %  - the final time
cada1f1 = length(tau);
N.f = cada1f1 - 1;
%User Line: N = length(tau)-1;
cada1f1 = N.f + 1;
cada1f2 = nstates*cada1f1;
stateIndices.f = 1:cada1f2;
%User Line: stateIndices = 1:nstates*(N+1);
cada1f1 = N.f + 1;
cada1f2 = nstates*cada1f1;
cada1f3 = cada1f2 + 1;
cada1f4 = N.f + 1;
cada1f5 = nstates*cada1f4;
cada1f6 = ncontrols*N.f;
cada1f7 = cada1f5 + cada1f6;
controlIndices.f = cada1f3:cada1f7;
%User Line: controlIndices = (nstates*(N+1)+1):(nstates*(N+1)+ncontrols*N);
cada1f1 = length(controlIndices.f);
cada1f2 = controlIndices.f(cada1f1);
t0Index.f = cada1f2 + 1;
%User Line: t0Index = controlIndices(end)+1;
tfIndex.f = t0Index.f + 1;
%User Line: tfIndex = t0Index+1;
stateVector.dz0 = z.dz0(Gator1Data.Index1);
stateVector.f = z.f(stateIndices.f);
%User Line: stateVector = z(stateIndices);
controlVector.dz0 = z.dz0(Gator1Data.Index2);
controlVector.f = z.f(controlIndices.f);
%User Line: controlVector = z(controlIndices);
t0.dz0 = z.dz0(1002);
t0.f = z.f(t0Index.f);
%User Line: t0 = z(t0Index);
tf.dz0 = z.dz0(1003);
tf.f = z.f(tfIndex.f);
%User Line: tf = z(tfIndex);
%User Line: % Reshape the state and control parts of the NLP decision vector
%User Line: % to matrices of size (N+1) by nstates and N by ncontrols, respectively.
cada1f1 = N.f + 1;
state.dz0 = stateVector.dz0;
state.f = reshape(stateVector.f,cada1f1,nstates);
%User Line: state   = reshape(stateVector,N+1,nstates);
control.dz0 = controlVector.dz0;
control.f = reshape(controlVector.f,N.f,ncontrols);
%User Line: control = reshape(controlVector,N,ncontrols);
cada1f1 = size(state.f,1);
cada1f2 = cada1f1 - 1;
cada1f3 = 1:cada1f2;
stateLGR.dz0 = state.dz0(Gator1Data.Index3);
stateLGR.f = state.f(cada1f3,:);
%User Line: stateLGR = state(1:end-1,:);
%User Line: % Identify the different components of the state column-wise from stateLGR
x.dz0 = stateLGR.dz0(Gator1Data.Index4);
x.f = stateLGR.f(:,1);
%User Line: x = stateLGR(:,1);
u.dz0 = control.dz0; u.f = control.f;
%User Line: u = control;
%User Line: % form cost
cada1td1 = zeros(2,1);
cada1td1(2) = tf.dz0;
cada1td1(1) = cada1td1(1) + -t0.dz0;
cada1f1dz0 = cada1td1;
cada1f1 = tf.f - t0.f;
cada1f2dz0 = cada1f1dz0./2;
cada1f2 = cada1f1/2;
cada1f3dz0 = 0.5.*cada1f2dz0;
cada1f3 = cada1f2*0.5;
cada1td1 = x.f(:).*x.dz0;
cada1td1 = cada1td1 + x.f(:).*x.dz0;
cada1f4dz0 = cada1td1;
cada1f4 = x.f.*x.f;
cada1td1 = u.f(:).*u.dz0;
cada1td1 = cada1td1 + u.f(:).*u.dz0;
cada1f5dz0 = cada1td1;
cada1f5 = u.f.*u.f;
cada1td1 = zeros(1000,1);
cada1td1(Gator1Data.Index5) = cada1f4dz0;
cada1td1(Gator1Data.Index6) = cada1td1(Gator1Data.Index6) + cada1f5dz0;
cada1f6dz0 = cada1td1;
cada1f6 = cada1f4 + cada1f5;
cada1tempdz0 = cada1f3dz0(Gator1Data.Index7);
cada1tf1 = cada1f6(Gator1Data.Index8);
cada1td1 = zeros(2000,1);
cada1td1(Gator1Data.Index9) = cada1tf1(:).*cada1tempdz0;
cada1td1(Gator1Data.Index10) = cada1td1(Gator1Data.Index10) + cada1f3.*cada1f6dz0;
Lagrangian.dz0 = cada1td1;
Lagrangian.f = cada1f3*cada1f6;
%User Line: Lagrangian = ((tf-t0)/2)*0.5*(x.*x+u.*u);
cada1f1 = psStuff.w.';
cada1td1 = sparse(Gator1Data.Index11,Gator1Data.Index12,Lagrangian.dz0,500,1002);
cada1td1 = cada1f1*cada1td1;
cada1td1 = cada1td1(:);
J.dz0 = full(cada1td1(Gator1Data.Index13));
J.f = cada1f1*Lagrangian.f;
%User Line: J   = psStuff.w.'*Lagrangian;
obj.dz0 = J.dz0; obj.f = J.f;
%User Line: obj = J;
obj.dz0_size = 1003;
obj.dz0_location = Gator1Data.Index14;
end


function ADiGator_LoadData()
global ADiGator_hyperSensitiveObj_ADiGatorJac
ADiGator_hyperSensitiveObj_ADiGatorJac = load('hyperSensitiveObj_ADiGatorJac.mat');
return
end