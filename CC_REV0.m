%%
%Leitura dos dados de entrada
clc
clear
%BARRAS CONECTADDAS A UG ALTERAR NUMERO PARA ZERO
%ORDENAR DADOS DE LINHA PARA QUE p<q SEMPRE
%ORDENAR DADOS DA BARRA(:,1) DO MENOR PARA O MAIOR
%input formato [p q Zpq]
input=[0 1 0.15j]; %0 - 1 0.15j
input=[input;1 2 0.15j]; %0 - 1 0.15j
input=[input;2 3 0.15j]; %0 - 1 0.15j
input=[input;3 4 0.15j]; %0 - 1 0.15j
%%
%Definição limites da matriz Zbarra
barras=unique(input(:,[1 2]));
z_barra=zeros(size(barras,1)-1,size(barras,1)-1); %retirar barra de ref
%%
%Algorítmo de construção da matriz Zbarra
for m=1:1:size(input,1)
    if(input(m,1)==0) %conexão com barra ref   
        barra=input(m,2);
        if(z_barra(barra,barra)==0) %barra é uma nova conexão
            disp('barra ainda não adicionada')
            z_barra(barra,barra)=input(m,3);
        else
            disp('barra ja adicionada')
        end
    end
    if(input(m,1)~=0) %não conexão com barra ref
        p=input(m,1);
        q=input(m,2);
        if(z_barra(q,q)==0) %barra q é uma nova conexão
            disp('barra ainda não adicionada')
            z_barra(q,q)=z_barra(p,p)+input(m,3);
        else
            disp('barra ja adicionada')
        end
    end

end







