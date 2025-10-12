%% test state_preparation
clc ; clear ; close all ; 
% addpath('QCLAB');

logging = true; % no record 
circuit_sim = false ; 

disp('Testing f3(n) for n = 1 to 15:');
for n = 1:10
    N = pow2(n) ; 
    state_complex = randn(N,1) + randn(N,1) .* 1j ; 
    state_complex = state_complex ./ norm(state_complex,2) ; 
    [circuit, global_phase, CNOT_count] = state_preparation( state_complex, 1, logging ) ; 
    fprintf('fun_CNOT(%d) = %.4f       CNOT(%d) %.4f  \n', n,  number_CNOT_state_preparation(n), n, CNOT_count) ; 
end 



%% CNOT's Upper bound of PB method for state preparation

n = 15 ;

lower_bound2 = @(n) ceil( (1/2).*pow2(n)-(3/4).*n-(1/4) );
lower_bound2(n)

% CNOT's lower bound (even case) function derived from 
% Martin Plesch and ˇCaslav Brukner. Quantum-state preparation with universal gate decompositions.
% Phys. Rev. A, 83:032302, Mar 2011 
lower_bound = @(n) ceil((1/2).*pow2(n)-(1/2).*n-(1/2)) ;
lower_bound(n)

% CNOT's upper bound (even case) function derived from Section V. A
% Martin Plesch and ˇCaslav Brukner. Quantum-state preparation with universal gate decompositions.
% Phys. Rev. A, 83:032302, Mar 2011 
upper_bound_PB_even = @(n) pow2(n./2) - 1 + (23/24).*pow2(n) - (3/2).*pow2(n./2+1) + 8/3 ; 
upper_bound_PB_even(n)

% CNOT's upper bound (odd case) function derived from Section V. B
% Martin Plesch and ˇCaslav Brukner. Quantum-state preparation with universal gate decompositions.
% Phys. Rev. A, 83:032302, Mar 2011 
upper_bound_PB_odd = @(n) pow2((n-1)./2) - 1 + (23/48).*pow2((n-1)) - (3/2).*pow2((n-1)./2) + 4/3 + (23/48).*pow2((n-1)+2) - (3/2).*pow2((n-1)./2+1) + 4/3 ;
upper_bound_PB_odd(n)

%% State Praprartion with complex amplitude

logging = true ; 

for n = 2:9
    N = pow2(n) ;
    for k = 1:5
        state_complex = randn(N,1) + randn(N,1) .* 1j ; 
        state_complex = state_complex ./ norm(state_complex,2) ; 
        [circuit, global_phase, CNOT_count] = state_preparation( state_complex, 1, logging ) ; 
        % [circuit2, global_phase2, CNOT_count2] = state_preparation2( state_complex, 1, logging ) ; 
        % [circuit, global_phase] = state_preparation( state_complex, 0, true ) ; 

        M = circuit.matrix ;
        circuit_state = M(:,1) .* global_phase ; 
        fprintf(" %2d      %4d      %5s   \n" ,n,CNOT_count,norm(circuit_state - state_complex) ) ; 

        % M2 = circuit2.matrix ;
        % circuit_state2 = M2(:,1) .* global_phase2 ; 
        % fprintf(" %2d      %4d      %5s \n" ,n,CNOT_count2,norm(circuit_state2 - state_complex) ) ;

        if norm(circuit_state - state_complex) > 1e-10 
            circuit_state 
            state_complex 
            circuit_state - state_complex 
            break ; 
        end
    end
end


%% test quantum state preparation by Schmidt decomposition  
% Reference: 
%   1. M. Plesch and C. Brukner, Quantum-state preparation with 
%      universal gate decompositions, Phys. Rev. A 83, 032302 (2011).  

clc; clear; 
logging  = true ; 

n = 7 ; 
N = pow2(n) ; 
state_complex = randn(N,1) + randn(N,1) .* 1j ; 
state_complex = state_complex ./ norm(state_complex,2) ; 


circuit0 = qclab.QCircuit( n ); 

reshape_vec = reshape( state_complex, pow2(floor(n/2)), [] ) ; 
[U, S, V] = svd( reshape_vec.' ) ; 

DU = diag( exp( randn(pow2(ceil(n/2)),1).*1j ) ) ; 
DV = diag( exp( randn(pow2(floor(n/2)),1).*1j ) ) ; 


[circuitS, global_phase, ~] = state_preparation( diag(DU'*S*conj(DV)), ceil(n/2) - floor(n/2), logging ) ;
circuit0.push_back(circuitS) ;   

M = circuit0.matrix .* global_phase ; 
state_simulation = zeros( pow2(n), 1 ) ; 
matrix_eye = eye( pow2( ceil(n/2) ) ) ; 

for i = 0 : floor(n/2) - 1 
    circuit0.push_back( qclab.qgates.CNOT(i + ceil(n/2) - floor(n/2), i + ceil(n/2)  ) ) ; 
end 

is_last_U = true ; 
circuit1 = qclab.QCircuit( ceil(n/2), 0 ) ; 
[circuit1, global_phase0, ~, ~] = recursive_decouple_quantum_circuit(circuit1, U*DU, ceil(n/2), 1, eye(4), is_last_U, true ) ; 
global_phase = global_phase .* global_phase0 ; 

circuit2 = qclab.QCircuit( floor(n/2), ceil(n/2) ) ; 
[circuit2, global_phase0, ~, ~] = recursive_decouple_quantum_circuit(circuit2, conj(V)*DV, floor(n/2), 1, eye(4), is_last_U, true ) ; 
global_phase = global_phase .* global_phase0 ; 

circuit0.push_back( circuit1 ) ; 
circuit0.push_back( circuit2 ) ; 

circuit0.draw()
M = circuit0.matrix .* global_phase ;
norm( M(:,1) - state_complex ) 

%% test unitary synthesis 

n = 4 ; 
circuit = qclab.QCircuit( n ); 
U = generate_U2n( n ) ; 
[circuit, global_phase, Delta, CNOT_count] = recursive_decouple_quantum_circuit( circuit, U, n, 1, eye(4), false, true ) ; 
norm(  kron(eye(pow2(n-2)), Delta') * circuit.matrix .* global_phase - U ) 

circuit = qclab.QCircuit( n ); 
[circuit, global_phase, Delta, CNOT_count] = recursive_decouple_quantum_circuit(circuit, U, n, 1, eye(4), true, true ) ; 
norm(  circuit.matrix .* global_phase - U ) 

%% test quantum circuit of isometry


n = 5 ; 
U = generate_U2n( ceil(n/2) ) ; 
circuit = qclab.QCircuit( ceil(n/2) ); 
[circuit, global_phase, Delta, CNOT_count] = isometry( U', circuit, ceil(n/2), true ) ;  
C = circuit.ctranspose ; 
M = C.matrix * kron(eye(pow2(ceil(n/2)-2)), Delta) .* global_phase' ; 
norm( M(:,1:pow2(floor(n/2))) - U(:,1:pow2(floor(n/2))) ) 


function [circuit, global_phase, Delta, CNOT_count] = isometry(U, circuit, n, logging) 
% Prepare an pow2(floor(n/2)) * pow2(floor(n/2)) unitary whose first pow2(floor(n/2)) rows are the isometry form pow2(n) to pow2(floor(n/2)) 
%   but attention that above n is different from the input of isometry() 
% Input: 
%   U            -- unitary to be decoupled, whose whose first pow2(floor(n/2)) columns are the isometry form pow2(ceil(n/2)) to pow2(floor(n/2)) 
%   n            -- parameter 
%   Delta        -- a diagonal matrix size of 4*4 in the bottom of a quantum ciruict 
%   is_last_U    -- boolean logical that whether U is the last unitary in the decomposition 
% --------------------------------------------------- 
% Output: 
%   circuit      -- qclab circuit 
%   global_phase -- the global phase of the circuit 
%   Delta        -- a diagonal matrix size of 4*4 in the bottom of a quantum ciruict 
%   CNOT_count   -- the CNOT count 
% ---------------------------------------------------
% test program: 
    % n = 5 ; 
    % U = generate_U2n( ceil(n/2) ) ; 
    % circuit = qclab.QCircuit( ceil(n/2) ); 
    % [circuit, global_phase, Delta, CNOT_count] = isometry( U', circuit, ceil(n/2), true ) ;  
    % C = circuit.ctranspose ; 
    % M = C.matrix * kron(eye(pow2(ceil(n/2)-2)), Delta) .* global_phase' ; 
    % norm( M(:,1:pow2(floor(n/2))) - U(:,1:pow2(floor(n/2))) ) 
% --------------------------------------------------- 

    assert( logging == true ) ; 
    if logging, CNOT_count = 0 ; end
    
    % Step 1. Compute the rotation parameters 
    [~, B1, B2, C] = transfomed_ZXZ(U) ; 
    % [A, B1, B2, C] = transfomed_ZXZ(U) ; 
    % U = [eye(pow2(n-1)), zeros(pow2(n-1)); zeros(pow2(n-1)), A]*kron(H,eye(pow2(n-1))) * [B1,zeros(pow2(n-1));zeros(pow2(n-1)),B2] * kron(H, eye(pow2(n-1))) * [ eye(pow2(n-1)), zeros(pow2(n-1));zeros(pow2(n-1)), C ] ;
    [WC, DC, VC] = decouple_unitary( eye(pow2(n-1)), C ) ; 
    % norm([WC, zeros(pow2(n-1)); zeros(pow2(n-1)), WC] * [DC, zeros(pow2(n-1)); zeros(pow2(n-1)), DC'] * [VC, zeros(pow2(n-1)); zeros(pow2(n-1)), VC] - [eye(pow2(n-1)), zeros(pow2(n-1)); zeros(pow2(n-1)), C] )  
    Z = [1 0; 0 -1] ; 
    [WB, DB, VB] = decouple_unitary( B1*WC, B2*WC*kron(Z,eye(pow2(n-2))) ) ; 
    % norm([WB, zeros(pow2(n-1)); zeros(pow2(n-1)), WL] * [DB, zeros(pow2(n-1)); zeros(pow2(n-1)), DB'] * [VB, zeros(pow2(n-1)); zeros(pow2(n-1)), VB] - [B1*WC, zeros(pow2(n-1)); zeros(pow2(n-1)), B2*WC*kron(Z,eye(pow2(n-2)))] )  
    
    % Step 2. Construct quantum circuit 
    global_phase = 1 ; Delta = eye(4) ; 
    [circuit, global_phase, Delta, CNOT_count0] = recursive_decouple_quantum_circuit(circuit, VC, n-1, global_phase, Delta, false, logging) ; 
    if logging, CNOT_count = CNOT_count0 + CNOT_count ; end
    % M = kron(eye(2), VC) ;
    % norm( kron( eye(pow2(n-2)), Delta') * circuit.matrix .* global_phase -  M )  
    [circuit, CNOT_count0] = uniformly_controlled_rotation_z(circuit, DC, n, 'L', logging) ; 
    if logging, CNOT_count = CNOT_count0 + CNOT_count ; end 
    % M = [DC, zeros(pow2(n-1)); zeros(pow2(n-1)), DC'] * M ; 
    % norm( kron( eye(pow2(n-2)), Delta') * kron(CNOT21, eye(pow2(n-2)))  * circuit.matrix .* global_phase -  M )  
    circuit.push_back( qclab.qgates.Hadamard(circuit.nbQubits - n) ) ; 
    [circuit, global_phase, Delta, CNOT_count0] = recursive_decouple_quantum_circuit(circuit, VB, n-1, global_phase, Delta, false, logging) ; 
    if logging, CNOT_count = CNOT_count0 + CNOT_count ; end
    [circuit, CNOT_count0] = uniformly_controlled_rotation_z(circuit, DB, n, 'M', logging) ; 
    if logging, CNOT_count = CNOT_count0 + CNOT_count ; end
    [circuit, global_phase, Delta, CNOT_count0] = recursive_decouple_quantum_circuit(circuit, WB, n-1, global_phase, Delta, false, logging) ; 
    if logging, CNOT_count = CNOT_count0 + CNOT_count ; end
    circuit.push_back( qclab.qgates.Hadamard(circuit.nbQubits - n) ) ;

end 


function [U] = generate_U2n(n)
    [U, ~] = qr( randn(pow2(n),pow2(n)) + 1i*randn(pow2(n),pow2(n)) );
    % U = Q / (det(Q)^(1/pow2(n))); % From U(2^n) to SU(2^n)
end


