function [A,B,C,D] = valato(index,I_asta,K_asta,C_smorz)
% matrice delle masse o delle inerzie
%  valori esempio per provare se funge
% index = 7;
% I_asta = 3;
% K_asta = 1000;
% C_smorz = 0.05;
for i = 1:index % indice in input
m(i,i) = I_asta; % kgm2


    % K matrice delle rigidezze
    % stesso indice dato di simulink
    if i == 1
        k(i,i)      = K_asta;
        k(i,i+1)    = -K_asta;
        c(i,i)      = C_smorz;
        c(i,i+1)    = -C_smorz;
    elseif i ~= Naste %dato in ingresso il numero di aste
        k(i,i-1)    = -K_asta;
        k(i,i)      = 2*K_asta;
        k(i,i+1)    = -K_asta;

        c(i,i-1)    = -C_smorz;
        c(i,i)      = 2*C_smorz;
        c(i,i+1)    = -C_smorz;
    else
        k(i,i-1)    = -K_asta;
        k(i,i)      = K_asta;

        c(i,i-1)    = -C_smorz;
        c(i,i)      = C_smorz;
    end
end

minv = inv(m);
cinv = -minv*c;
kinv = -minv*k;
C = eye(2*index);
D = [zeros(index);zeros(index)];

if m*minv == eye(index)
    disp('tutt''apposto')
else
    error('checcazzaccio hai fatto?')
end

for j = 1:index*2
    for i = 1:index*2

        if i <= index && j <= index
            A(i,j) = 0;
            B(i,j) = 0;
        elseif i<= index && j> index
            if j == i+index
                A(i,j) = 1;
            else
                A(i,j) = 0;
            end
        elseif i > index && j <= index
            A(i,j) = kinv(i-index,j);
            B(i,j) = minv(i-index,j);
        elseif i > index && j > index
            A(i,j) = cinv(i-index,j-index);
        end
    end
end

end