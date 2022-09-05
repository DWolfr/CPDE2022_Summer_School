% Exercise for the Summer School 
% A Practical Introduction to Control, Numerics and Machine Learning
% on the IFAC CPDE 2022 
% Workshop on Control of Systems Governed by Partial Differential Equations
% Dr. Daniel Veldman (d.w.m.veldman@gmail.com)

function T = OCP_compute_temperature(time, T0, E, A, Bd, B, U, LU)

% the last input argument is an LU factorization of E-dt/2*A, which you can
% use to speed up computations

T = zeros(length(T0), length(time));
T(:,1) = T0;
for ii = 1:length(time)-1
    dt = time(ii+1) - time(ii);

    % Security for ukm1
    if ii == 1
        ukm1 = [0;0];
    else
        ukm1 = U(:,ii);
    end

    % Calculate RHS
    rhs = ((E + (dt/2) * A) * T(:,ii)) +  dt * B * ((U(:,ii) - ukm1)/2) + dt * Bd;
    
    if nargin >= 8
        % LU factorization is available. 
        T(:,ii+1) = LU \ rhs;
    else
        % LU factorization is not available. (START WITH THIS LINE)
        T(:,ii+1) = (E - (dt/2) * A)\rhs;
    end
end