%%
%Leitura dos dados de entrada
clc
clear
%BARRAS CONECTADDAS A UG ALTERAR NUMERO PARA ZERO
%ENTRAR COM OS SEGUINTES DADOS MANUALMENTE:
rf=0; %RESISTÊNCIA DE FALTA [pu]
Sb=100; %DEFINIR POTENCIA DE BASE DO SISTEMA [MVA]
barra_cc=4; %NOME DA BARRA COMO NO ARQUIVO DE ENTRADA
%AS COMPONENTES DE SEQUÊNCIAS ESTÃO REFERIDAS A FASE A
%CURTO CIRCUITO BIFÁSICO: CURTO NAS FASES B E C
%CURTO CIRCUITO MONOFASICO: FASE A À TERRA
%%
%ordem das 11 colunas:
%p|q|z+[pu ou ohm/km]|z0[pu ou ohm/km]|L[km]|Vpb[kV]|Vqb[kV]|Sb_equip[MVA]|Vp_equip[kV]|Vq_equip[kV]|tipo_seq_0
%%%%%%%%%%%%DADOS DE ENTRADA%%%%%%%%%%%%%%%%%%%
entrada=csvread('./4bar_radial_OK.txt'); %ex livro cap 6 cc barra 4
% entrada=csvread('./5bar_OK.txt'); %ex livro cap 6 cc barra 5
% entrada=csvread('./8bar_ex14_OK.txt'); %ex 14 lista p1, não há seq zero nos dados monofásico inválido cc barra 7
% entrada=csvread('./8bar_ex30_OK.txt'); %ex 30 lista p1 cc barra 7
% entrada=csvread('./4bar_p1ex1_OK.txt'); %p1 ex1 matriz zbarra
% entrada=csvread('./8bar_p1ex2_OK.txt'); %p1 ex2 matriz zbarra
% entrada=csvread('./8bar_p1ex2_sem_ramo_6_4_OK.txt'); %p1 ex2 itens a b c cc barra 6
% entrada=csvread('./IEEE13.txt'); %cc barra 
%03/11/19:todos os resultados dos arquivos testados acima batem com os resultados de referencia 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
barraNome=unique(entrada(:,[1 2])); %esta ordenado menor para o maior devido função unique
barraNum=0:1:length(barraNome)-1; %barra 0 é sempre o primeiro elemento do vetor
V_barra_base=[entrada(:,1) entrada(:,6);entrada(:,2) entrada(:,7)]; %tensão de base por barra, podem haver registros repetidos
barraNum_sem_barras_ficticias=barraNum; %usado para calculo de tensão nas barras
%ORDENAR DADOS DE LINHA PARA QUE p<q
entrada(:,[1 2])=sort(entrada(:,[1 2]),2);
%ORDENAR DADOS DA BARRA(:,1) DO MENOR PARA O MAIOR
entrada=sortrows(entrada,1);
for i=1:1:length(V_barra_base)
    V_barra_base(i,1)=barraNum(find(barraNome==V_barra_base(i,1))); %mapear nome barra ao numero
end
for i=1:1:length(entrada(:,1))
    entrada(i,1)=barraNum(find(barraNome==entrada(i,1))); %mapear nome barra p ao numero
    entrada(i,2)=barraNum(find(barraNome==entrada(i,2))); %mapear nome barra q ao numero
end
entrada_sem_ramos_ficticios=entrada; %usado para calculo de corrente nos ramos
entrada_aux=[];
for m=1:1:size(entrada,1) %tratamento de topologia seq pos e zero
    X=entrada(m,:); %para facilitar leitura das variaveis
    [p,q,zp,z0,L,Vpb,Vqb,Sequip,Vp_equip,Vq_equip,t_z0]=deal(X(1),X(2),X(3),X(4),X(5),X(6),X(7),X(8),X(9),X(10),X(11));
    if(L>0) %LT
        zp=zp*L*(Sb/Vpb^2);
        z0=z0*L*(Sb/Vpb^2);
        entrada_aux=[entrada_aux;[p q zp z0 L Vpb Vqb Sb Vp_equip Vq_equip t_z0]];        
    elseif(Sequip>0) %TR e UGs
        if(Vp_equip==0) %evitar casos ramo p ou q é ref para geradores
            zp=zp*((Vq_equip^2)/Sequip)*(Sb/Vqb^2);
            z0=zp;
        else
            zp=zp*((Vp_equip^2)/Sequip)*(Sb/Vpb^2);
            z0=zp;
        end 
        if(t_z0==1) %trafo seq zero conexão p-ficticia-0 |q
           barraNum=[barraNum max(barraNum)+1];
           ficticia=barraNum(end);
           entrada_aux=[entrada_aux;[p ficticia zp/2 z0/2 L Vpb Vqb Sb Vp_equip Vq_equip t_z0]]; %conexão p-ficticia
           entrada_aux=[entrada_aux;[0 ficticia 9999 z0/2 L Vpb Vqb Sb Vp_equip Vq_equip t_z0]]; %conexão ficticia-0
           entrada_aux=[entrada_aux;[ficticia q zp/2 9999 L Vpb Vqb Sb Vp_equip Vq_equip t_z0]]; %conexão ficticia-q        
        elseif(t_z0==2) %trafo seq zero conexão p| 0-ficticia-q 
           barraNum=[barraNum max(barraNum)+1];
           ficticia=barraNum(end);
           entrada_aux=[entrada_aux;[p ficticia zp/2 9999 L Vpb Vqb Sb Vp_equip Vq_equip t_z0]]; %conexão p-ficticia
           entrada_aux=[entrada_aux;[0 ficticia 9999 z0/2 L Vpb Vqb Sb Vp_equip Vq_equip t_z0]]; %conexão ficticia-0
           entrada_aux=[entrada_aux;[ficticia q zp/2 z0/2 L Vpb Vqb Sb Vp_equip Vq_equip t_z0]]; %conexão ficticia-q 
        elseif(t_z0==3) %trafo seq zero conexão direta p-q
           entrada_aux=[entrada_aux;[p q zp z0 L Vpb Vqb Sb Vp_equip Vq_equip t_z0]]; %conexão p-q
        elseif(t_z0==4) %trafo seq zero conexão p| |q 
           entrada_aux=[entrada_aux;[p q zp 9999 L Vpb Vqb Sb Vp_equip Vq_equip t_z0]]; %conexão p desconectado de p seq 0
        elseif(t_z0==5) %trafo seq zero conexão p-ficticia-q e ficticia-0 (nucleo envolvido Yaterrado-Yaterrado
           barraNum=[barraNum max(barraNum)+1];
           ficticia=barraNum(end);
           entrada_aux=[entrada_aux;[p ficticia zp/2 z0/2 L Vpb Vqb Sb Vp_equip Vq_equip t_z0]]; %conexão p-ficticia
           entrada_aux=[entrada_aux;[0 ficticia 9999 5*z0 L Vpb Vqb Sb Vp_equip Vq_equip t_z0]]; %conexão ficticia-0
           entrada_aux=[entrada_aux;[ficticia q zp/2 z0/2 L Vpb Vqb Sb Vp_equip Vq_equip t_z0]]; %conexão ficticia-q
        elseif(t_z0==6) %trafo seq zero conexão p-ficticia-q e ficticia-0 (nucleo envolvido Yaterrado-Yaterrado
           barraNum=[barraNum max(barraNum)+1];
           ficticia=barraNum(end);
           entrada_aux=[entrada_aux;[p ficticia zp/2 0.4*z0 L Vpb Vqb Sb Vp_equip Vq_equip t_z0]]; %conexão p-ficticia
           entrada_aux=[entrada_aux;[0 ficticia 9999 0.45*z0 L Vpb Vqb Sb Vp_equip Vq_equip t_z0]]; %conexão ficticia-0
           entrada_aux=[entrada_aux;[ficticia q zp/2 9999 L Vpb Vqb Sb Vp_equip Vq_equip t_z0]]; %conexão ficticia-q         
        elseif(t_z0==7) %trafo seq zero conexão p-ficticia-q e ficticia-0 (nucleo envolvido Yaterrado-Yaterrado
           barraNum=[barraNum max(barraNum)+1];
           ficticia=barraNum(end);
           entrada_aux=[entrada_aux;[p ficticia zp/2 9999 L Vpb Vqb Sb Vp_equip Vq_equip t_z0]]; %conexão p-ficticia
           entrada_aux=[entrada_aux;[0 ficticia 9999 0.45*z0 L Vpb Vqb Sb Vp_equip Vq_equip t_z0]]; %conexão ficticia-0
           entrada_aux=[entrada_aux;[ficticia q zp/2 0.4*z0 L Vpb Vqb Sb Vp_equip Vq_equip t_z0]]; %conexão ficticia-q
        end
    else %equivalente de sistemam em pu
        entrada_aux=[entrada_aux;[p q zp z0 L Vpb Vqb Sb Vp_equip Vq_equip t_z0]];
    end
end
entrada=entrada_aux;
%%
%Definição limites da matriz Zbarra
total_barras_sem_ficticias=length(barraNum_sem_barras_ficticias)-1; %retirar barra de ref
total_barras=length(barraNum)-1; %retirar barra de ref
z_pos_barra=zeros(total_barras,total_barras);
z_zero_barra=zeros(total_barras,total_barras);
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
            z_pos_barra=z_pos_barra-(1/(z_pos_barra(q,q)+z_pos_lt))*z_pos_barra(:,q)*z_pos_barra(q,:);
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
            z_zero_barra=z_zero_barra-(1/(z_zero_barra(q,q)+z_zero_lt))*z_zero_barra(:,q)*z_zero_barra(q,:);
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
disp(['====RESULTADO PARA CURTO-CIRCUITO BARRA ', num2str(barra_cc), '====='])
barra_cc=barraNum(find(barraNome==barra_cc));
%%
%TRIFÁSICO
%Cálculo corrente de curto-circuito
[row,~]=find(V_barra_base(:,1)==barra_cc);
Vbase=V_barra_base(row(1),2);
Ibase=Sb/(Vbase*sqrt(3));
Icc3f=Ibase/z_pos_barra(barra_cc,barra_cc);
Icc3f_abs=abs(Icc3f);
Icc3f_angle=radtodeg(angle(Icc3f));
%Cálculo tensão nas barras
V3f=zeros(1,total_barras_sem_ficticias);
for i=1:1:total_barras_sem_ficticias %quantidade de barras menos a referencia e barras ficticias
    [row,~]=find(V_barra_base(:,1)==i);
    Vbase=V_barra_base(row(1),2)/sqrt(3);
    V3f(i)=Vbase*(1-z_pos_barra(i,barra_cc)/z_pos_barra(barra_cc,barra_cc));
    V3f_abs(i)=abs(V3f(i));
    V3f_angle(i)=radtodeg(angle(V3f(i)));
end
%Cálculo fluxo de corrente nos ramos
Ipq3f=zeros(1,size(entrada_sem_ramos_ficticios,1));
for m=1:1:size(entrada_sem_ramos_ficticios,1)
    X=entrada_sem_ramos_ficticios(m,:); %para facilitar leitura das variaveis
    [p,q,zp,z0,L,Vpb,Vqb,Sequip,Vp_equip,Vq_equip,t_z0]=deal(X(1),X(2),X(3),X(4),X(5),X(6),X(7),X(8),X(9),X(10),X(11));
    if(L>0) %LT
        zp=zp*L*(Sb/Vpb^2);
    elseif(Sequip>0) %TR e UGs
        if(Vp_equip==0) %evitar casos ramo p ou q é ref para geradores
            zp=zp*((Vq_equip^2)/Sequip)*(Sb/Vqb^2);
        else
            zp=zp*((Vp_equip^2)/Sequip)*(Sb/Vpb^2);
        end 
    end
    if(p==0)
        Ibase=Sb/(Vqb*sqrt(3)); %sempre considerando corrente do ramo p ao ramo q
        Ipq3f(m)=Ibase*z_pos_barra(q,barra_cc)/(z_pos_barra(barra_cc,barra_cc)*zp);
        Ipq3f_abs(m)=abs(Ipq3f(m));
        Ipq3f_angle(m)=radtodeg(angle(Ipq3f(m)));
    else
        Ibase=Sb/(Vpb*sqrt(3)); %sempre considerando corrente do ramo p ao ramo q
        Ipq3f(m)=Ibase*(z_pos_barra(q,barra_cc)-z_pos_barra(p,barra_cc))/(z_pos_barra(barra_cc,barra_cc)*zp);
        Ipq3f_abs(m)=abs(Ipq3f(m));
        Ipq3f_angle(m)=radtodeg(angle(Ipq3f(m)));
    end  
end
%%
%BIFÁSICO
%Cálculo corrente de curto-circuito
[row,~]=find(V_barra_base(:,1)==barra_cc);
Vbase=V_barra_base(row(1),2);
Ibase=Sb/(Vbase*sqrt(3));
Icc2fA=0;
Icc2fA_abs=abs(Icc2fA);
Icc2fA_angle=radtodeg(angle(Icc2fA));
Icc2fB=Ibase*-1j*sqrt(3)*(1/(2*z_pos_barra(barra_cc,barra_cc)));
Icc2fB_abs=abs(Icc2fB);
Icc2fB_angle=radtodeg(angle(Icc2fB));
Icc2fC=-Icc2fB;
Icc2fC_abs=abs(Icc2fC);
Icc2fC_angle=radtodeg(angle(Icc2fC));
%Cálculo tensão nas barras
for i=1:1:total_barras_sem_ficticias %quantidade de barras menos a referencia e barras ficticias    
    [row,~]=find(V_barra_base(:,1)==i);
    Vbase=V_barra_base(row(1),2)/sqrt(3);
    V2fA(i)=Vbase*1;
    V2fA_abs(i)=abs(V2fA(i));
    V2fA_angle(i)=radtodeg(angle(V2fA(i)));
    V2fB(i)=Vbase*(exp(4j*pi/3)+1j*sqrt(3)*z_pos_barra(i,barra_cc)/(2*z_pos_barra(barra_cc,barra_cc)));
    V2fB_abs(i)=abs(V2fB(i));
    V2fB_angle(i)=radtodeg(angle(V2fB(i)));
    V2fC(i)=Vbase*(exp(2j*pi/3)-1j*sqrt(3)*z_pos_barra(i,barra_cc)/(2*z_pos_barra(barra_cc,barra_cc)));
    V2fC_abs(i)=abs(V2fC(i));
    V2fC_angle(i)=radtodeg(angle(V2fC(i)));
end
%Cálculo fluxo de corrente nos ramos
for m=1:1:size(entrada_sem_ramos_ficticios,1)
    X=entrada_sem_ramos_ficticios(m,:); %para facilitar leitura das variaveis
    [p,q,zp,z0,L,Vpb,Vqb,Sequip,Vp_equip,Vq_equip,t_z0]=deal(X(1),X(2),X(3),X(4),X(5),X(6),X(7),X(8),X(9),X(10),X(11));
    if(L>0) %LT
        zp=zp*L*(Sb/Vpb^2);
    elseif(Sequip>0) %TR e UGs
        if(Vp_equip==0) %evitar casos ramo p ou q é ref para geradores
            zp=zp*((Vq_equip^2)/Sequip)*(Sb/Vqb^2);
        else
            zp=zp*((Vp_equip^2)/Sequip)*(Sb/Vpb^2);
        end 
    end
    if(p==0)
        Ibase=Sb/(Vqb*sqrt(3)); %sempre considerando corrente do ramo p ao ramo q
        Ipq2fA(m)=0;
        Ipq2fA_abs(m)=abs(Ipq2fA(m));
        Ipq2fA_angle(m)=radtodeg(angle(Ipq2fA(m)));
        Ipq2fB(m)=Ibase*1j*sqrt(3)*(-z_pos_barra(q,barra_cc))/(2*z_pos_barra(barra_cc,barra_cc)*zp);
        Ipq2fB_abs(m)=abs(Ipq2fB(m));
        Ipq2fB_angle(m)=radtodeg(angle(Ipq2fB(m)));
        Ipq2fC(m)=Ibase*1j*sqrt(3)*(z_pos_barra(q,barra_cc))/(2*z_pos_barra(barra_cc,barra_cc)*zp);
        Ipq2fC_abs(m)=abs(Ipq2fC(m));
        Ipq2fC_angle(m)=radtodeg(angle(Ipq2fC(m)));
    else
        Ibase=Sb/(Vpb*sqrt(3)); %sempre considerando corrente do ramo p ao ramo q
        Ipq2fA(m)=0;
        Ipq2fA_abs(m)=abs(Ipq2fA(m));
        Ipq2fA_angle(m)=radtodeg(angle(Ipq2fA(m)));
        Ipq2fB(m)=Ibase*1j*sqrt(3)*(z_pos_barra(p,barra_cc)-z_pos_barra(q,barra_cc))/(2*z_pos_barra(barra_cc,barra_cc)*zp);
        Ipq2fB_abs(m)=abs(Ipq2fB(m));
        Ipq2fB_angle(m)=radtodeg(angle(Ipq2fB(m)));
        Ipq2fC(m)=Ibase*1j*sqrt(3)*(z_pos_barra(q,barra_cc)-z_pos_barra(p,barra_cc))/(2*z_pos_barra(barra_cc,barra_cc)*zp);
        Ipq2fC_abs(m)=abs(Ipq2fC(m));
        Ipq2fC_angle(m)=radtodeg(angle(Ipq2fC(m)));
    end
end
%%
%MONOFÁSICO
%Cálculo corrente de curto-circuito
[row,~]=find(V_barra_base(:,1)==barra_cc);
Vbase=V_barra_base(row(1),2);
Ibase=Sb/(Vbase*sqrt(3));
Icc1fA=Ibase*3/(2*z_pos_barra(barra_cc,barra_cc)+z_zero_barra(barra_cc,barra_cc)+3*rf);
Icc1fA_abs=abs(Icc1fA);
Icc1fA_angle=radtodeg(angle(Icc1fA));
Icc1fB=0;
Icc1fB_abs=abs(Icc1fB);
Icc1fB_angle=radtodeg(angle(Icc1fB));
Icc1fC=0;
Icc1fC_abs=abs(Icc1fC);
Icc1fC_angle=radtodeg(angle(Icc1fC));
%Cálculo tensão nas barras
for i=1:1:total_barras_sem_ficticias %quantidade de barras menos a referencia e barras ficticias    
    [row,~]=find(V_barra_base(:,1)==i);
    Vbase=V_barra_base(row(1),2)/sqrt(3);
    V1fA(i)=Vbase*(1-((2*z_pos_barra(i,barra_cc)+z_zero_barra(i,barra_cc))/(2*z_pos_barra(barra_cc,barra_cc)+z_zero_barra(barra_cc,barra_cc)+3*rf)));
    V1fA_abs(i)=abs(V1fA(i));
    V1fA_angle(i)=radtodeg(angle(V1fA(i)));
    V1fB(i)=Vbase*(exp(4j*pi/3)-((z_zero_barra(i,barra_cc)-z_pos_barra(i,barra_cc))/(2*z_pos_barra(barra_cc,barra_cc)+z_zero_barra(barra_cc,barra_cc)+3*rf)));
    V1fB_abs(i)=abs(V1fB(i));
    V1fB_angle(i)=radtodeg(angle(V1fB(i)));
    V1fC(i)=Vbase*(exp(2j*pi/3)-((z_zero_barra(i,barra_cc)-z_pos_barra(i,barra_cc))/(2*z_pos_barra(barra_cc,barra_cc)+z_zero_barra(barra_cc,barra_cc)+3*rf)));
    V1fC_abs(i)=abs(V1fC(i));
    V1fC_angle(i)=radtodeg(angle(V1fC(i)));
end
%Cálculo fluxo de corrente nos ramos
for m=1:1:size(entrada_sem_ramos_ficticios,1)
    X=entrada_sem_ramos_ficticios(m,:); %para facilitar leitura das variaveis
    [p,q,zp,z0,L,Vpb,Vqb,Sequip,Vp_equip,Vq_equip,t_z0]=deal(X(1),X(2),X(3),X(4),X(5),X(6),X(7),X(8),X(9),X(10),X(11));
    if(L>0) %LT
        zp=zp*L*(Sb/Vpb^2);
        z0=z0*L*(Sb/Vpb^2);
    elseif(Sequip>0) %TR e UGs
        if(Vp_equip==0) %evitar casos ramo p ou q é ref para geradores
            zp=zp*((Vq_equip^2)/Sequip)*(Sb/Vqb^2);
            z0=z0*((Vq_equip^2)/Sequip)*(Sb/Vqb^2);
        else
            zp=zp*((Vp_equip^2)/Sequip)*(Sb/Vpb^2);
            z0=z0*((Vp_equip^2)/Sequip)*(Sb/Vpb^2);
        end 
    end
    if(p==0)
        Ibase=Sb/(Vqb*sqrt(3)); %sempre considerando corrente do ramo p ao ramo q
        Ipq1fA(m)=Ibase*(1/(2*z_pos_barra(barra_cc,barra_cc)+z_zero_barra(barra_cc,barra_cc)+3*rf))*(2*((z_pos_barra(q,barra_cc))/zp)+((z_zero_barra(q,barra_cc))/z0));
        Ipq1fA_abs(m)=abs(Ipq1fA(m));
        Ipq1fA_angle(m)=radtodeg(angle(Ipq1fA(m)));
        Ipq1fB(m)=Ibase*(1/(2*z_pos_barra(barra_cc,barra_cc)+z_zero_barra(barra_cc,barra_cc)+3*rf))*(((-z_pos_barra(q,barra_cc))/zp)+((z_zero_barra(q,barra_cc))/z0));
        Ipq1fB_abs(m)=abs(Ipq1fB(m));
        Ipq1fB_angle(m)=radtodeg(angle(Ipq1fB(m)));
        Ipq1fC(m)=Ipq1fB(m);
        Ipq1fC_abs(m)=abs(Ipq1fC(m));
        Ipq1fC_angle(m)=radtodeg(angle(Ipq1fC(m)));
    else 
        Ibase=Sb/(Vpb*sqrt(3)); %sempre considerando corrente do ramo p ao ramo q
        Ipq1fA(m)=Ibase*(1/(2*z_pos_barra(barra_cc,barra_cc)+z_zero_barra(barra_cc,barra_cc)+3*rf))*(2*((z_pos_barra(q,barra_cc)-z_pos_barra(p,barra_cc))/zp)+((z_zero_barra(q,barra_cc)-z_zero_barra(p,barra_cc))/z0));
        Ipq1fA_abs(m)=abs(Ipq1fA(m));
        Ipq1fA_angle(m)=radtodeg(angle(Ipq1fA(m)));
        Ipq1fB(m)=Ibase*(1/(2*z_pos_barra(barra_cc,barra_cc)+z_zero_barra(barra_cc,barra_cc)+3*rf))*(((z_pos_barra(p,barra_cc)-z_pos_barra(q,barra_cc))/zp)+((z_zero_barra(q,barra_cc)-z_zero_barra(p,barra_cc))/z0));
        Ipq1fB_abs(m)=abs(Ipq1fB(m));
        Ipq1fB_angle(m)=radtodeg(angle(Ipq1fB(m)));
        Ipq1fC(m)=Ipq1fB(m);
        Ipq1fC_abs(m)=abs(Ipq1fC(m));
        Ipq1fC_angle(m)=radtodeg(angle(Ipq1fC(m)));
    end
end

%%
%MAPEAMENTO NOME DAS BARRAS COMO NOS DADOS DE ENTRADA
for i=1:1:length(entrada_sem_ramos_ficticios(:,1))
    entrada_sem_ramos_ficticios(i,1)=barraNome(find(barraNum==entrada_sem_ramos_ficticios(i,1))); %mapear numero barra p ao nome
    entrada_sem_ramos_ficticios(i,2)=barraNome(find(barraNum==entrada_sem_ramos_ficticios(i,2))); %mapear numero barra q ao nome
end
%SAÍDA DE DADOS TRIFÁSICO
disp(['*****CURTO-CIRCUITO TRIFÁSICO*****'])
TI3F=array2table([barraNome(find(barraNum==barra_cc)) Icc3f_abs Icc3f_angle], 'VariableNames',{'Barra', 'I_cc_kA', 'deg'});
disp(TI3F)
TV3F=array2table([barraNome(2:end) V3f_abs.' V3f_angle.'], 'VariableNames',{'Barra', 'Vfase_kV', 'deg'});
disp(TV3F)
TIPQ3F=array2table([entrada_sem_ramos_ficticios(:,[1 2]) Ipq3f_abs.' Ipq3f_angle.'], 'VariableNames',{'p', 'q', 'I_kA', 'deg'});
disp(TIPQ3F)
writetable(TI3F,'I_cc_3F.txt','Delimiter',' ');
writetable(TV3F,'V_cc_3F.txt','Delimiter',' ');
writetable(TIPQ3F,'Ipq_cc_3F.txt','Delimiter',' ');
disp(['***************************************'])
%SAÍDA DE DADOS BIFÁSICO
disp(['*****CURTO-CIRCUITO BIFÁSICO******'])
TI2F=array2table([barraNome(find(barraNum==barra_cc)) Icc2fA_abs Icc2fA_angle Icc2fB_abs Icc2fB_angle Icc2fC_abs Icc2fC_angle], 'VariableNames',{'Barra', 'IfA_cc_kA', 'deg_fA', 'IfB_cc_kA', 'deg_fB', 'IfC_cc_kA', 'deg_fC'});
disp(TI2F)
TV2F=array2table([barraNome(2:end) V2fA_abs.' V2fA_angle.' V2fB_abs.' V2fB_angle.' V2fC_abs.' V2fC_angle.'], 'VariableNames',{'Barra', 'VfA_kV', 'VfA_deg', 'VfB_kV', 'VfB_deg', 'VfC_kV', 'VfC_deg'});
disp(TV2F)
TIPQ2F=array2table([entrada_sem_ramos_ficticios(:,[1 2]) Ipq2fA_abs.' Ipq2fA_angle.' Ipq2fB_abs.' Ipq2fB_angle.' Ipq2fC_abs.' Ipq2fC_angle.'], 'VariableNames',{'p', 'q', 'IfA_kA', 'deg_fA', 'IfB_kA', 'deg_fB', 'IfC_kA', 'deg_fC'});
disp(TIPQ2F)
writetable(TI2F,'I_cc_2F.txt','Delimiter',' ');
writetable(TV2F,'V_cc_2F.txt','Delimiter',' ');
writetable(TIPQ2F,'Ipq_cc_2F.txt','Delimiter',' ');
disp(['***************************************'])
%SAÍDA DE DADOS MONOFÁSICO
disp(['*****CURTO-CIRCUITO MONOFÁSICO******'])
TI1F=array2table([barraNome(find(barraNum==barra_cc)) Icc1fA_abs Icc1fA_angle Icc1fB_abs Icc1fB_angle Icc1fC_abs Icc1fC_angle], 'VariableNames',{'Barra', 'IfA_cc_kA', 'deg_fA', 'IfB_cc_kA', 'deg_fB', 'IfC_cc_kA', 'deg_fC'});
disp(TI1F)
TV1F=array2table([barraNome(2:end) V1fA_abs.' V1fA_angle.' V1fB_abs.' V1fB_angle.' V1fC_abs.' V1fC_angle.'], 'VariableNames',{'Barra', 'VfA_kV', 'VfA_deg', 'VfB_kV', 'VfB_deg', 'VfC_kV', 'VfC_deg'});
disp(TV1F)
TIPQ1F=array2table([entrada_sem_ramos_ficticios(:,[1 2]) Ipq1fA_abs.' Ipq1fA_angle.' Ipq1fB_abs.' Ipq1fB_angle.' Ipq1fC_abs.' Ipq1fC_angle.'], 'VariableNames',{'p', 'q', 'IfA_kA', 'deg_fA', 'IfB_kA', 'deg_fB', 'IfC_kA', 'deg_fC'});
disp(TIPQ1F)
writetable(TI1F,'I_cc_1F.txt','Delimiter',' ');
writetable(TV1F,'V_cc_1F.txt','Delimiter',' ');
writetable(TIPQ1F,'Ipq_cc_1F.txt','Delimiter',' ');
disp(['***************************************'])
