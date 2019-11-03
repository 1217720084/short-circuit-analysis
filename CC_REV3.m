%%
%Leitura dos dados de entrada
clc
clear
%BARRAS CONECTADDAS A UG ALTERAR NUMERO PARA ZERO
%TOPOLOGIA z0 CONECTAR TERRA A BARRA 0
%%
%ordem colunas: p   | q     | z+    | z0
stream=fopen('./5bar.txt');
data=textscan(stream, '%s', 'Delimiter', ' ');
fclose(stream);
p={data{1}{1:4:length(data{1})}};
q={data{1}{2:4:length(data{1})}};
z_pos={data{1}{3:4:length(data{1})}};
z_zero={data{1}{4:4:length(data{1})}};
barName=unique([p q]); %esta ordenado menor para o maior devido unique
barNum=[0 [1:1:length(barName)-1]]; %barra 0 é definida nos dados de 
for i=1:1:length(p)
    for j=1:1:length(barName)
        if strcmp(p(i),barName(j))
            p_num(i)=barNum(j); %adaptar p para algoritmo
        end
    end    
end
for i=1:1:length(q)
    for j=1:1:length(barName)
        if strcmp(q(i),barName(j))
            q_num(i)=barNum(j); %adaptar q para algoritmo
        end
    end    
end
for i=1:1:length(z_pos)
   z_pos_num(i)=str2num(z_pos{i}); %convesão de str para num
end
for i=1:1:length(z_zero)
   z_zero_num(i)=str2num(z_zero{i}); %convesão de str para num
end
input=[p_num.' q_num.' z_pos_num.' z_zero_num.'];
%ORDENAR DADOS DE LINHA PARA QUE p<q
input(:,[1 2])=sort(input(:,[1 2]),2);
%ORDENAR DADOS DA BARRA(:,1) DO MENOR PARA O MAIOR
input=sortrows(input,1);
%%
%Definição limites da matriz Zbarra
barras=unique(input(:,[1 2]));
barras_validas=size(barras,1)-1; %retirar barra de ref
z_pos_barra=zeros(barras_validas,barras_validas);
z_zero_barra=zeros(barras_validas,barras_validas);
%%
%Algorítmo de construção da matriz Zbarra de sequencia positiva 
for m=1:1:size(input,1)
    p=input(m,1);
    q=input(m,2);
    z_pos_lt=input(m,3);
    disp(['conexão p-', num2str(p), ' q-', num2str(q), ':'])
    if(p==0) %conexão com barra ref   
        if(z_pos_barra(q,q)==0) %barra é uma nova conexão
            disp(['*barra ', num2str(q), ' ainda não adicionada'])
            z_pos_barra(q,q)=z_pos_lt;
        else
            disp(['*barra ', num2str(q), ' já adicionada'])           
            z_pos_barra=z_pos_barra-(1/z_pos_barra(q,q)+z_pos_lt)*z_pos_barra(:,q)*z_pos_barra(q,:);
        end
    else
        if(z_pos_barra(q,q)==0) %barra q é uma nova conexão
            disp(['*barra ', num2str(q), ' ainda não adicionada'])
            z_pos_barra(:,q)=z_pos_barra(:,p);
            z_pos_barra(q,:)=z_pos_barra(p,:);
            z_pos_barra(q,q)=z_pos_barra(p,p)+z_pos_lt;
        else
            disp(['*barra ', num2str(q), ' já adicionada'])
            z_pos_barra=z_pos_barra-(1/(z_pos_barra(q,q)+z_pos_barra(p,p)-(2*z_pos_barra(p,q))+z_pos_lt))*(z_pos_barra(:,p)-z_pos_barra(:,q))*(z_pos_barra(p,:)-z_pos_barra(q,:));
        end
    end
end
%%
%Algorítmo de construção da matriz Zbarra de sequencia zero
for m=1:1:size(input,1)
    p=input(m,1);
    q=input(m,2);
    z_zero_lt=input(m,4);
    disp(['conexão p-', num2str(p), ' q-', num2str(q), ':'])
    if(p==0) %conexão com barra ref   
        if(z_zero_barra(q,q)==0) %barra é uma nova conexão
            disp(['*barra ', num2str(q), ' ainda não adicionada'])
            z_zero_barra(q,q)=z_zero_lt;
        else
            disp(['*barra ', num2str(q), ' já adicionada'])           
            z_zero_barra=z_zero_barra-(1/z_zero_barra(q,q)+z_zero_lt)*z_zero_barra(:,q)*z_zero_barra(q,:);
        end
    else
        if(z_zero_barra(q,q)==0) %barra q é uma nova conexão
            disp(['*barra ', num2str(q), ' ainda não adicionada'])
            z_zero_barra(:,q)=z_zero_barra(:,p);
            z_zero_barra(q,:)=z_zero_barra(p,:);
            z_zero_barra(q,q)=z_zero_barra(p,p)+z_zero_lt;
        else
            disp(['*barra ', num2str(q), ' já adicionada'])
            z_zero_barra=z_zero_barra-(1/(z_zero_barra(q,q)+z_zero_barra(p,p)-(2*z_zero_barra(p,q))+z_zero_lt))*(z_zero_barra(:,p)-z_zero_barra(:,q))*(z_zero_barra(p,:)-z_zero_barra(q,:));
        end
    end
end
%%
%ESCOLHER BARRA EM CC
barra_cc=5;
%%
%TRIFÁSICO
%Cálculo corrente de curto-circuito
Icc=1/z_pos_barra(barra_cc,barra_cc);
array2table([barra_cc Icc], 'VariableNames',{'Barra', 'I_cc'})
%Cálculo tensão nas barras
V=zeros(1,barras_validas);
for i=1:1:barras_validas %quantidade de barras menos a referencia
    V(i)=1-z_pos_barra(i,barra_cc)/z_pos_barra(barra_cc,barra_cc);
end
array2table([[1:1:barras_validas];V].', 'VariableNames',{'Barra', 'V_pu'})
%Cálculo fluxo de corrente nos ramos
[row,col]=find(input(:,1)==0);
ramos_exceto_referencia=input(row(end)+1:1:size(input,1),:);
I=zeros(1,size(ramos_exceto_referencia,1));
for m=1:1:size(ramos_exceto_referencia,1)
    p=ramos_exceto_referencia(m,1);
    q=ramos_exceto_referencia(m,2);
    z_pos_lt=ramos_exceto_referencia(m,3);
    I(m)=(z_pos_barra(q,barra_cc)-z_pos_barra(p,barra_cc))/(z_pos_barra(barra_cc,barra_cc)*z_pos_lt);
end
array2table([ramos_exceto_referencia(:,[1 2]) I.'], 'VariableNames',{'p', 'q', 'I_pu'})

