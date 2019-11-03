%%
%Leitura dos dados de entrada
clc
clear
%BARRAS CONECTADDAS A UG ALTERAR NUMERO PARA ZERO
%TOPOLOGIA z0 CONECTAR TERRA A BARRA 0
%%
%ordem colunas: p   | q     | z+    | z0
entrada=csvread('./8bar.txt');
barraNome=unique(entrada(:,[1 2])); %esta ordenado menor para o maior devido função unique
barraNum=0:1:length(barraNome)-1; %barra 0 é sempre o primeiro elemento do vetor
for i=1:1:length(entrada(:,1))
    entrada(i,1)=barraNum(find(barraNome==entrada(i,1))); %mapear nome barra p ao numero
end
for i=1:1:length(entrada(:,2))
    entrada(i,2)=barraNum(find(barraNome==entrada(i,2))); %mapear nome barra q ao numero
end
%ORDENAR DADOS DE LINHA PARA QUE p<q
entrada(:,[1 2])=sort(entrada(:,[1 2]),2);
%ORDENAR DADOS DA BARRA(:,1) DO MENOR PARA O MAIOR
entrada=sortrows(entrada,1);
for m=1:1:size(entrada,1)
end
%%
%Definição limites da matriz Zbarra
barras_validas=size(barraNome,1)-1; %retirar barra de ref
z_pos_barra=zeros(barras_validas,barras_validas);
z_zero_barra=zeros(barras_validas,barras_validas);
%%
%Algorítmo de construção da matriz Zbarra de sequencia positiva 
for m=1:1:size(entrada,1)
    p=entrada(m,1);
    q=entrada(m,2);
    z_pos_lt=entrada(m,3);
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
        elseif(z_pos_barra(p,p)==0) %barra p é uma nova conexão
            disp(['*barra ', num2str(p), ' ainda não adicionada'])
            z_pos_barra(:,p)=z_pos_barra(:,q);
            z_pos_barra(p,:)=z_pos_barra(q,:);
            z_pos_barra(p,p)=z_pos_barra(q,q)+z_pos_lt;
        else
            disp(['*barra ', num2str(q), ' já adicionada'])
            z_pos_barra=z_pos_barra-(1/(z_pos_barra(q,q)+z_pos_barra(p,p)-(2*z_pos_barra(p,q))+z_pos_lt))*(z_pos_barra(:,p)-z_pos_barra(:,q))*(z_pos_barra(p,:)-z_pos_barra(q,:));
        end
    end
end
%%
%Algorítmo de construção da matriz Zbarra de sequencia zero
for m=1:1:size(entrada,1)
    p=entrada(m,1);
    q=entrada(m,2);
    z_zero_lt=entrada(m,4);
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
        elseif(z_zero_barra(p,p)==0) %barra p é uma nova conexão
            disp(['*barra ', num2str(p), ' ainda não adicionada'])
            z_zero_barra(:,p)=z_zero_barra(:,q);
            z_zero_barra(p,:)=z_zero_barra(q,:);
            z_zero_barra(p,p)=z_zero_barra(q,q)+z_zero_lt;
        else
            disp(['*barra ', num2str(q), ' já adicionada'])
            z_zero_barra=z_zero_barra-(1/(z_zero_barra(q,q)+z_zero_barra(p,p)-(2*z_zero_barra(p,q))+z_zero_lt))*(z_zero_barra(:,p)-z_zero_barra(:,q))*(z_zero_barra(p,:)-z_zero_barra(q,:));
        end
    end
end
%%
%ESCOLHER BARRA EM CC
barra_cc=605; %nome como no arquivo de entrada
disp(['====RESULTADO PARA CURTO-CIRCUITO BARRA ', num2str(barra_cc), '====='])
barra_cc=barraNum(find(barraNome==barra_cc));
%%
%TRIFÁSICO
%Cálculo corrente de curto-circuito
Icc3f=1/z_pos_barra(barra_cc,barra_cc);
%Cálculo tensão nas barras
V3f=zeros(1,barras_validas);
for i=1:1:barras_validas %quantidade de barras menos a referencia
    V3f(i)=1-z_pos_barra(i,barra_cc)/z_pos_barra(barra_cc,barra_cc);
end
%Cálculo fluxo de corrente nos ramos
[row,col]=find(entrada(:,1)==0);
ramos_exceto_referencia=entrada(row(end)+1:1:size(entrada,1),:); %ignorar ramos ligados a barra 0
Ipq3f=zeros(1,size(ramos_exceto_referencia,1));
for m=1:1:size(ramos_exceto_referencia,1)
    p=ramos_exceto_referencia(m,1);
    q=ramos_exceto_referencia(m,2);
    z_pos_lt=ramos_exceto_referencia(m,3);
    Ipq3f(m)=(z_pos_barra(q,barra_cc)-z_pos_barra(p,barra_cc))/(z_pos_barra(barra_cc,barra_cc)*z_pos_lt);
end
%%
%Saída de dados
disp(array2table([barraNome(find(barraNum==barra_cc)) Icc3f], 'VariableNames',{'Barra', 'I_cc'}))
disp(array2table([barraNome(2:end) V3f.'], 'VariableNames',{'Barra', 'V_pu'}))
for i=1:1:length(entrada(:,1))
    entrada(i,1)=barraNome(find(barraNum==entrada(i,1))); %mapear numero barra p ao nome
end
for i=1:1:length(entrada(:,2))
    entrada(i,2)=barraNome(find(barraNum==entrada(i,2))); %mapear numero barra q ao nome
end
ramos_exceto_referencia=entrada(row(end)+1:1:size(entrada,1),:); %ignorar ramos ligados a barra 0
disp(array2table([ramos_exceto_referencia(:,[1 2]) Ipq3f.'], 'VariableNames',{'p', 'q', 'I_pu'}))

