function [MPCparams] = ParametersMPCMulti(A,B,C,Q,R,iu,PTv)

% Determine size of matrices
nState = size(A,1);         % number of states

% number of tracked outputs
nOut = size(C,1);

BB = B(:,iu);             % extract B matrix for corresponding control input
nCont = length(iu);         % number of control input variations
D  = zeros(nOut,nCont);     % initialise direct control effect matrix D

nPH = (size(PTv,2) - 1);

F = zeros(nPH*nOut,nState);
G = zeros(nPH*nOut,nCont*(nPH)); % exclude D;

for nn = 1:nPH

    for pp = 1:nOut

        F(pp + (nn-1) * nOut,:) = C(pp,:)*expm(A*PTv(pp,nn+1));	
		
    end

end

for nn = 1:nPH      % building row wise of G

    for pp = 1:nOut     % for each output

        for mm = 1:nn       % building colomn wise of G            
		
            G(pp + (nn-1) * nOut,[1+(mm-1)*nCont:nCont+(mm-1)*nCont]) = C(pp,:) * expm(A*(PTv(pp,nn+1)-PTv(pp,mm+1))) * ...
                                                                        pinv(A) * (expm(A*(PTv(pp,mm+1)-PTv(pp,mm)))-eye(size(A)))*BB;

        end

    end

end

Gt = G';                        % transpose G
M = inv(Gt*Q*G + R);          % element from cost function minimisation

H = Gt*Q*G + R;
f = -Gt*Q;

% formulate augmented stability matrix

K = M * Gt * Q;

MPCparams.F = F;  % free response matrix
MPCparams.G = G;  % forced response matrix
MPCparams.M = M;
MPCparams.K = K;  % optimised gain
MPCparams.Q = Q;  % tracking error weights
MPCparams.R = R;  % control weights
MPCparams.nPH = nPH;  % number of coincidence points
MPCparams.nOut = nOut;    % number of outputs
MPCparams.A = A;  % state sensitivity matrix
MPCparams.B = BB; % control sensitivity matrix
MPCparams.C = C;
MPCparams.PTv = PTv;  % coincidence points
% MPCparams.dT = dT;
MPCparams.QR = norm(Q)/norm(R); % Q to R ratio (Q/R  > 1e-5 good, Q/R < 1e-5,-6 bad)

MPCparams.H = H;  % hessian
MPCparams.f = f;

return

function E = expmeig(A)

[V,D] = eig(A);
E = V * diag(exp(diag(D))) / V;

E = real(E);

return