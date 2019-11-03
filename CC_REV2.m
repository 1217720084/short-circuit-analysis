%%
%Leitura dos dados de entrada
clc
clear
%BARRAS CONECTADDAS A UG ALTERAR NUMERO PARA ZERO
%ORDENAR DADOS DE LINHA PARA QUE p<q SEMPRE
%ORDENAR DADOS DA BARRA(:,1) DO MENOR PARA O MAIOR
%input formato [p q Zpq]
%%
%teste sistema radial
% input=[0 1 0.15j]; %0 - 1 0.15j
% input=[input;1 2 0.15j]; %0 - 1 0.15j
% input=[input;2 3 0.15j]; %0 - 1 0.15j
% input=[input;3 4 0.15j]; %0 - 1 0.15j
%%
%teste sistema exercicio lista
input=[0 1 0.2j]; %0 - 1 0.2j
input=[input;0 2 0.3j]; %0 - 2 0.3j
input=[input;1 2 0.8j]; %1 - 2 0.8j
input=[input;1 3 0.8j]; %1 - 3 0.8j
input=[input;2 4 1j]; %2 - 4 1j
input=[input;3 4 0.6j]; %3 - 4 0.8j
%%
%Definição limites da matriz Zbarra
barras=unique(input(:,[1 2]));
barras_validas=size(barras,1)-1; %retirar barra de ref
z_barra=zeros(barras_validas,barras_validas);
%%
%Algorítmo de construção da matriz Zbarra
for m=1:1:size(input,1)
    p=input(m,1);
    q=input(m,2);
    zlt=input(m,3);
    disp(['conexão p-', num2str(p), ' q-', num2str(q), ':'])
    if(p==0) %conexão com barra ref   
        if(z_barra(q,q)==0) %barra é uma nova conexão
            disp(['*barra ', num2str(q), ' ainda não adicionada'])
            z_barra(q,q)=zlt;
        else
            disp(['*barra ', num2str(q), ' já adicionada'])           
            z_barra=z_barra-(1/z_barra(q,q)+zlt)*z_barra(:,q)*z_barra(q,:);
        end
    else
        if(z_barra(q,q)==0) %barra q é uma nova conexão
            disp(['*barra ', num2str(q), ' ainda não adicionada'])
            z_barra(:,q)=z_barra(:,p);
            z_barra(q,:)=z_barra(p,:);
            z_barra(q,q)=z_barra(p,p)+zlt;
        else
            disp(['*barra ', num2str(q), ' já adicionada'])
            z_barra=z_barra-(1/(z_barra(q,q)+z_barra(p,p)-(2*z_barra(p,q))+zlt))*(z_barra(:,p)-z_barra(:,q))*(z_barra(p,:)-z_barra(q,:));
        end
    end
end
%%
%ESCOLHER BARRA EM CC
barra_cc=2; 
%%
%Cálculo corrente de curto-circuito
Icc=1/z_barra(barra_cc,barra_cc);
array2table([barra_cc Icc], 'VariableNames',{'Barra', 'I_cc'})
%%
%Cálculo tensão nas barras
V=zeros(1,barras_validas);
for i=1:1:barras_validas %quantidade de barras menos a referencia
    V(i)=1-z_barra(i,barra_cc)/z_barra(barra_cc,barra_cc);
end
array2table([[1:1:barras_validas];V].', 'VariableNames',{'Barra', 'V_pu'})
%%
%Cálculo fluxo de corrente nos ramos
[row,col]=find(input(:,1)==0);
ramos_exceto_referencia=input(row(end)+1:1:size(input,1),:);
I=zeros(1,size(ramos_exceto_referencia,1));
for m=1:1:size(ramos_exceto_referencia,1)
    p=ramos_exceto_referencia(m,1);
    q=ramos_exceto_referencia(m,2);
    zlt=ramos_exceto_referencia(m,3);
    I(m)=(z_barra(q,barra_cc)-z_barra(p,barra_cc))/(z_barra(barra_cc,barra_cc)*zlt);
end
array2table([ramos_exceto_referencia(:,[1 2]) I.'], 'VariableNames',{'p', 'q', 'I_pu'})

