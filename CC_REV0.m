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
%Defini��o limites da matriz Zbarra
barras=unique(input(:,[1 2]));
z_barra=zeros(size(barras,1)-1,size(barras,1)-1); %retirar barra de ref
%%
%Algor�tmo de constru��o da matriz Zbarra
for m=1:1:size(input,1)
    if(input(m,1)==0) %conex�o com barra ref   
        barra=input(m,2);
        if(z_barra(barra,barra)==0) %barra � uma nova conex�o
            disp('barra ainda n�o adicionada')
            z_barra(barra,barra)=input(m,3);
        else
            disp('barra ja adicionada')
        end
    end
    if(input(m,1)~=0) %n�o conex�o com barra ref
        p=input(m,1);
        q=input(m,2);
        if(z_barra(q,q)==0) %barra q � uma nova conex�o
            disp('barra ainda n�o adicionada')
            z_barra(q,q)=z_barra(p,p)+input(m,3);
        else
            disp('barra ja adicionada')
        end
    end

end







