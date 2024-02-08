clear
clc
M=8;
CDE2S=3*10^2;
sigma=1.1;
r1=45;
RD=1000;

df1=1.5;
df2=2.25;

Theta1=0.12;
Theta2=0.10;
beta=0.002;
Lambda=[0.7,1.2,1.7,2.2];
Alpha1=-1;
Alpha2=20;
Alpha3=60;
Alpha4=-80;
c=50;
tD = zeros(1,200);
L = length(tD);
for i=1:L
    tD(i)=0.001*10^(9*i/99);
end

v = zeros(1,M);
for j=1:M
    
    v(j)=0;
    for k=fix((j+1)/2):1:min(j,M/2)
        v(j)=v(j) + k^(M/2) * factorial(2*k+1) / (factorial(M/2-k+1) * factorial(k+1) * factorial(k) * factorial(j-k+1) * factorial(2*k-j+1));%¸Ä½ø
        
    end
    v(j)=(-1)^(M/2+j)*v(j);
end

for k=1:4
    for i=1:L
        ft(i,k)=0;
        u=log(2)/tD(i);
        for j=1:M
            z=j*u;
            
            v1=1-(df1/(2+Theta1));
            v2=1-(df2/(2+Theta2));
            
            A=(r1^(1-(df1-Theta1)/2))*posimn(v1,v1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S));
            B=(r1^(-(df1-Theta1)/2))*((2+Theta1-df1)*posimn(v1,v1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S))+sqrt(z/CDE2S)*(r1^((2+Theta1)/2))*posimn(v1,v1+1,1,r1^((Theta1+2)/2),...
                (2/(Theta1+2))*sqrt(z/CDE2S)));
            C=(r1^(1-(df1-Theta1)/2))*((2+Theta1-df1)*posimn(v1,v1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S))-sqrt(z/CDE2S)*(1^((2+Theta1)/2))*posimn(v1+1,v1,1,r1^((Theta1+2)/2),...
                (2/(Theta1+2))*sqrt(z/CDE2S)));
            D=(r1^(-(df1-Theta1)/2))*(((2+Theta1-df1)^2)*posimn(v1,v1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S))-sqrt(z/CDE2S)*(2+Theta1-df1)*posimn(v1+1,v1,1,r1^((Theta1+2)/2),...
                (2/(Theta1+2))*sqrt(z/CDE2S))+sqrt(z/CDE2S)*(2+Theta1-df1)*(r1^(1+Theta1/2))*posimn(v1,v1+1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S))-(z/CDE2S)*(r1^(1+Theta1/2))...
                *posimn(v1+1,v1+1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S)));
            
            E=((r1*RD)^(1-(df2-Theta2)/2))*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S));
            F=(r1^(1-(df2-Theta2)/2))*(RD^(-(df2-Theta2)/2))*((2+Theta2-df2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))+sqrt(sigma*z/CDE2S)*(RD^(1+Theta2/2))...
                *posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
            G=(r1^(-(df2-Theta2)/2))*(RD^(1-(df2-Theta2)/2))*((2+Theta2-df2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-sqrt(sigma*z/CDE2S)*(r1^(1+Theta2/2))...
                *posimn(v2+1,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
            H=((r1*RD)^(-(df2-Theta2)/2))*(((2+Theta2-df2)^2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-sqrt(sigma*z/CDE2S)*(2+Theta2-df2)*(r1^(1+Theta2/2))...
                *posimn(v2+1,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))+sqrt(sigma*z/CDE2S)*(2+Theta2-df2)*(RD^(1+Theta2/2))*posimn(v2,v2+1,r1^(1+Theta2/2),...
                RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-(sigma*z/CDE2S)*((r1*RD)^(1+Theta2/2))*posimn(v2+1,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
            
            P=((r1*RD)^(1-(df2-Theta2)/2))*((v2/(2*z))*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))+(1/(2*sqrt(z)))*(2/(Theta2+2))*sqrt(sigma/CDE2S)...
                *(RD^(1+Theta2/2))*posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
            Q=(r1^(-(df2-Theta2)/2))*(RD^(1-(df2-Theta2)/2))*((v2/(2*z))*((2+Theta2-df2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-(r1^(1+Theta2/2))...
                *sqrt(sigma*z/CDE2S)*posimn(v2+1,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)))+(2/(Theta2+2))*sqrt(sigma/CDE2S)*(RD^(1+Theta2))*(1/(2*sqrt(z)))...
                *((2+Theta2-df2)*posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-(r1^(1+Theta2))*sqrt(sigma*z/CDE2S)*posimn(v2+1,v2+1,r1^(1+Theta2/2),...
                RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))));
  
            fz2=(Alpha3*(RD^2)+Alpha4*RD+c)*E-(Alpha1*RD+Alpha2)*P+RD*F;
            fm2=(Alpha3*(RD^2)+Alpha4*RD+c)*G-(Alpha1*RD+Alpha2)*Q+RD*H;
            fei2=fz2/fm2;
            fz1=Lambda(k)*r1^(df2-df1-(Theta2-Theta1))*A-fei2*B;
            fm1=Lambda(k)*r1^(df2-df1-(Theta2-Theta1))*C-fei2*D;
            fei1=fz1/fm1;
            ft(i,k) =  ft(i,k) + v(j) *(-fz1*beta/((z^2+z*beta)*fz1-z*fm1));
        end
        ft(i,k)=ft(i,k)*u;
        y(i,k)=-log(ft(i,k)+1)/beta;
    end
end

for k=1:4
    D_y(1,k) = tD(1) * (-3 * y(1,k) + 4 * y(2,k) - y(3,k)) / (4 * tD(2) - 3 * tD(1) - tD(3));
    D_y(2,k) = tD(2) * (y(3,k) - y(1,k)) / (tD(3) - tD(1));
    for i=3:L
        D_y(i,k) = tD(i) * (y(i - 2,k) - 4 * y(i - 1,k) + 3 * y(i,k)) / (tD(i - 2) - 4 * tD(i - 1) + 3 * tD(i));
    end
end

loglog(tD,y(:,1),'Color',[0.75, 0, 0.75], 'Linewidth',1.5)
hold on
loglog(tD,y(:,2),'Color',[0, 0.5, 0], 'Linewidth',1.5)
hold on
loglog(tD,y(:,3),'Color',[0.85, 0.33,0.2], 'Linewidth',1.5)
hold on
loglog(tD,y(:,4),'Color',[0, 0, 1], 'Linewidth',1.5)
hold on

loglog(tD,D_y(:,1),'--','Color',[0.75, 0, 0.75], 'Linewidth',1.5)
hold on
loglog(tD,D_y(:,2),'--','Color',[0, 0.5, 0], 'Linewidth',1.5)
hold on
loglog(tD,D_y(:,3),'--','Color',[0.85, 0.33,0.2], 'Linewidth',1.5)
hold on
loglog(tD,D_y(:,4),'--','Color',[0, 0, 1], 'Linewidth',1.5)
hold on

axis([10^(-1), 10^10, 10^(0), 10^3]);
hold on

hold on
xlabel('\itt_{D}','FontName','Times New Roman','FontSize',10) 
ylabel('\itP_{wD} , P\prime_{wD}' ,'FontName','Times New Roman','FontSize',10) 

h=legend('$$\lambda=0.7$$','$$\lambda=1.2$$','$$\lambda=1.7$$','$$\lambda=2.2$$');
set(h,'Interpreter','latex')

box off
grid on


