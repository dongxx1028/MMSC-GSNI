clear
clc
M=2;
CDE2S=3*10^2;
sigma=1.1;
r1=45;
RD=1000;
df1=1.5;
df2=2.25;
Theta1=0.12;
Theta2=0.10;
beta=0.002;
Lambda=1.2;
           
tD = zeros(1,200);
L = length(tD);
for i=1:L
    tD(i)=0.001*10^(9*i/99);
end

v = zeros(1,M);
for j=1:M
    
    v(j)=0;
    for k=fix((j+1)/2):1:min(j,M/2)
        v(j)=v(j) + k^(M/2) * factorial(2*k+1) / (factorial(M/2-k+1) * factorial(k+1) * factorial(k) * factorial(j-k+1) * factorial(2*k-j+1));
        
    end
    v(j)=(-1)^(M/2+j)*v(j);
end

for k=1:1
    Alpha1=0;
    Alpha2=0;
    Alpha3=0;
    Alpha4=0;
    c=0;
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
            F=(r1^(1-(df2-Theta2)/2))*(RD^(-(df2-Theta2)/2))*((2+Theta2-df2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))+sqrt(sigma*z/CDE2S)...
                *(RD^(1+Theta2/2))*posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
            G=(r1^(-(df2-Theta2)/2))*(RD^(1-(df2-Theta2)/2))*((2+Theta2-df2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-sqrt(sigma*z/CDE2S)...
                *(r1^(1+Theta2/2))*posimn(v2+1,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
            H=((r1*RD)^(-(df2-Theta2)/2))*(((2+Theta2-df2)^2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-sqrt(sigma*z/CDE2S)*(2+Theta2-df2)...
                *(r1^(1+Theta2/2))*posimn(v2+1,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))+sqrt(sigma*z/CDE2S)*(2+Theta2-df2)*(RD^(1+Theta2/2))...
                *posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-(sigma*z/CDE2S)*((r1*RD)^(1+Theta2/2))*posimn(v2+1,v2+1,r1^(1+Theta2/2),...
                RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
            
            P=((r1*RD)^(1-(df2-Theta2)/2))*((v2/(2*z))*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))+(1/(2*sqrt(z)))*(2/(Theta2+2))...
                *sqrt(sigma/CDE2S)*(RD^(1+Theta2/2))*posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
            Q=(r1^(-(df2-Theta2)/2))*(RD^(1-(df2-Theta2)/2))*((v2/(2*z))*((2+Theta2-df2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))...
                -(r1^(1+Theta2/2))*sqrt(sigma*z/CDE2S)*posimn(v2+1,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)))+(2/(Theta2+2))*sqrt(sigma/CDE2S)...
                *(RD^(1+Theta2))*(1/(2*sqrt(z)))*((2+Theta2-df2)*posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-(r1^(1+Theta2))...
                *sqrt(sigma*z/CDE2S)*posimn(v2+1,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))));

            fz2 = F;
            fm2 = H;
            fei2=fz2/fm2;
            
            fz1=Lambda*(r1^(df2-df1-(Theta2-Theta1)))*A-fei2*B;
            fm1=Lambda*(r1^(df2-df1-(Theta2-Theta1)))*C-fei2*D;
            fei1=fz1/fm1;
              
            ft(i,k) =  ft(i,k) + v(j) *(-fz1*beta(k)/((z^2+z*beta(k))*fz1-z*fm1));
        end
        ft(i,k)=ft(i,k)*u;
        y1(i,k)=-log(ft(i,k)+1)/beta(k);
    end
end

for k=1:1
    D_y1(1,k) = tD(1) * (-3 * y1(1,k) + 4 * y1(2,k) - y1(3,k)) / (4 * tD(2) - 3 * tD(1) - tD(3));
    D_y1(2,k) = tD(2) * (y1(3,k) - y1(1,k)) / (tD(3) - tD(1));
    for i=3:L
        D_y1(i,k) = tD(i) * (y1(i - 2,k) - 4 * y1(i - 1,k) + 3 * y1(i,k)) / (tD(i - 2) - 4 * tD(i - 1) + 3 * tD(i));
    end
end

for k=1:1
    Alpha1=0;
    Alpha2=0;
    Alpha3=0;
    Alpha4=0;
    c=50;
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
            F=(r1^(1-(df2-Theta2)/2))*(RD^(-(df2-Theta2)/2))*((2+Theta2-df2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))+sqrt(sigma*z/CDE2S)...
                *(RD^(1+Theta2/2))*posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
            G=(r1^(-(df2-Theta2)/2))*(RD^(1-(df2-Theta2)/2))*((2+Theta2-df2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-sqrt(sigma*z/CDE2S)...
                *(r1^(1+Theta2/2))*posimn(v2+1,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
            H=((r1*RD)^(-(df2-Theta2)/2))*(((2+Theta2-df2)^2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-sqrt(sigma*z/CDE2S)*(2+Theta2-df2)...
                *(r1^(1+Theta2/2))*posimn(v2+1,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))+sqrt(sigma*z/CDE2S)*(2+Theta2-df2)*(RD^(1+Theta2/2))...
                *posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-(sigma*z/CDE2S)*((r1*RD)^(1+Theta2/2))*posimn(v2+1,v2+1,r1^(1+Theta2/2),...
                RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
            
            P=((r1*RD)^(1-(df2-Theta2)/2))*((v2/(2*z))*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))+(1/(2*sqrt(z)))*(2/(Theta2+2))*sqrt(sigma/CDE2S)...
                *(RD^(1+Theta2/2))*posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
            Q=(r1^(-(df2-Theta2)/2))*(RD^(1-(df2-Theta2)/2))*((v2/(2*z))*((2+Theta2-df2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-(r1^(1+Theta2/2))...
                *sqrt(sigma*z/CDE2S)*posimn(v2+1,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)))+(2/(Theta2+2))*sqrt(sigma/CDE2S)*(RD^(1+Theta2))*(1/(2*sqrt(z)))...
                *((2+Theta2-df2)*posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-(r1^(1+Theta2))*sqrt(sigma*z/CDE2S)*posimn(v2+1,v2+1,r1^(1+Theta2/2),...
                RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))));

            fz2=E;
            fm2=G;
            fei2=fz2/fm2;
            
            fz1=Lambda*(r1^(df2-df1-(Theta2-Theta1)))*A-fei2*B;
            fm1=Lambda*(r1^(df2-df1-(Theta2-Theta1)))*C-fei2*D;
            fei1=fz1/fm1;
            
            
            ft(i,k) =  ft(i,k) + v(j) *(-fz1*beta(k)/((z^2+z*beta(k))*fz1-z*fm1));
        end
        ft(i,k)=ft(i,k)*u;
        y2(i,k)=-log(ft(i,k)+1)/beta(k);
    end
end


for k=1:1
    D_y2(1,k) = tD(1) * (-3 * y2(1,k) + 4 * y2(2,k) - y2(3,k)) / (4 * tD(2) - 3 * tD(1) - tD(3));
    D_y2(2,k) = tD(2) * (y2(3,k) - y2(1,k)) / (tD(3) - tD(1));
    for i=3:L
        D_y2(i,k) = tD(i) * (y2(i - 2,k) - 4 * y2(i - 1,k) + 3 * y2(i,k)) / (tD(i - 2) - 4 * tD(i - 1) + 3 * tD(i));
    end
end


for k=1:1
    
    for i=1:L
        ft(i,k)=0;
        u=log(2)/tD(i);
        for j=1:M
            z=j*u;
            
            v1=1-(df1/(2+Theta1));
            v2=1-(df2/(2+Theta2));
            
            A=(r1^(1-(df1-Theta1)/2))*posimn(v1,v1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S));
            B=(r1^(-(df1-Theta1)/2))*((2+Theta1-df1)*posimn(v1,v1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S))+sqrt(z/CDE2S)*(r1^((2+Theta1)/2))*posimn(v1,v1+1,1,r1^((Theta1+2)/2),(2/(Theta1+2))...
                *sqrt(z/CDE2S)));
            C=(r1^(1-(df1-Theta1)/2))*((2+Theta1-df1)*posimn(v1,v1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S))-sqrt(z/CDE2S)*(1^((2+Theta1)/2))*posimn(v1+1,v1,1,r1^((Theta1+2)/2),(2/(Theta1+2))...
                *sqrt(z/CDE2S)));
            D=(r1^(-(df1-Theta1)/2))*(((2+Theta1-df1)^2)*posimn(v1,v1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S))-sqrt(z/CDE2S)*(2+Theta1-df1)*posimn(v1+1,v1,1,r1^((Theta1+2)/2),(2/(Theta1+2))...
                *sqrt(z/CDE2S))+sqrt(z/CDE2S)*(2+Theta1-df1)*(r1^(1+Theta1/2))*posimn(v1,v1+1,1,r1^((Theta1+2)/2),(2/(Theta1+2))*sqrt(z/CDE2S))-(z/CDE2S)*(r1^(1+Theta1/2))...
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


            fz2=(r1^(1-(df2-Theta2)/2))*besselk(v2,(2/(Theta2+2))*sqrt(sigma*z/CDE2S)*(r1^(1+Theta2/2)));
            fm2=r1^(-((df2-Theta2)/2))*((2+Theta2-df2)*besselk(v2,(2/(Theta2+2))*sqrt(sigma*z/CDE2S)*(r1^(1+Theta2/2)))-sqrt(sigma*z/CDE2S)*(r1^(1+Theta2/2))*besselk(v2+1,(2/(Theta2+2))...
                *sqrt(sigma*z/CDE2S)*(r1^(1+Theta2/2))));
            fei2=fz2/fm2;
            
            fz1=Lambda*(r1^(df2-df1-(Theta2-Theta1)))*A-fei2*B;
            fm1=Lambda*(r1^(df2-df1-(Theta2-Theta1)))*C-fei2*D;
            fei1=fz1/fm1;
            
            
            ft(i,k) =  ft(i,k) + v(j) *(-fz1*beta(k)/((z^2+z*beta(k))*fz1-z*fm1));
        end
        ft(i,k)=ft(i,k)*u;
        y3(i,k)=-log(ft(i,k)+1)/beta(k);
    end
end


for k=1:1
    D_y3(1,k) = tD(1) * (-3 * y3(1,k) + 4 * y3(2,k) - y3(3,k)) / (4 * tD(2) - 3 * tD(1) - tD(3));
    D_y3(2,k) = tD(2) * (y3(3,k) - y3(1,k)) / (tD(3) - tD(1));
    for i=3:L
        D_y3(i,k) = tD(i) * (y3(i - 2,k) - 4 * y3(i - 1,k) + 3 * y3(i,k)) / (tD(i - 2) - 4 * tD(i - 1) + 3 * tD(i));
    end
end


for k=1:1
    Alpha1=-1;
    Alpha2=20;
    Alpha3=60;
    Alpha4=-80;
    c=50;
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
            F=(r1^(1-(df2-Theta2)/2))*(RD^(-(df2-Theta2)/2))*((2+Theta2-df2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))+sqrt(sigma*z/CDE2S)...
                *(RD^(1+Theta2/2))*posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
            G=(r1^(-(df2-Theta2)/2))*(RD^(1-(df2-Theta2)/2))*((2+Theta2-df2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-sqrt(sigma*z/CDE2S)...
                *(r1^(1+Theta2/2))*posimn(v2+1,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
            H=((r1*RD)^(-(df2-Theta2)/2))*(((2+Theta2-df2)^2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-sqrt(sigma*z/CDE2S)*(2+Theta2-df2)...
                *(r1^(1+Theta2/2))*posimn(v2+1,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))+sqrt(sigma*z/CDE2S)*(2+Theta2-df2)*(RD^(1+Theta2/2))...
                *posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-(sigma*z/CDE2S)*((r1*RD)^(1+Theta2/2))*posimn(v2+1,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),...
                (2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
            
            P=((r1*RD)^(1-(df2-Theta2)/2))*((v2/(2*z))*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))+(1/(2*sqrt(z)))*(2/(Theta2+2))*sqrt(sigma/CDE2S)...
                *(RD^(1+Theta2/2))*posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)));
            Q=(r1^(-(df2-Theta2)/2))*(RD^(1-(df2-Theta2)/2))*((v2/(2*z))*((2+Theta2-df2)*posimn(v2,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-(r1^(1+Theta2/2))...
                *sqrt(sigma*z/CDE2S)*posimn(v2+1,v2,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S)))+(2/(Theta2+2))*sqrt(sigma/CDE2S)*(RD^(1+Theta2))*(1/(2*sqrt(z)))...
                *((2+Theta2-df2)*posimn(v2,v2+1,r1^(1+Theta2/2),RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))-(r1^(1+Theta2))*sqrt(sigma*z/CDE2S)*posimn(v2+1,v2+1,r1^(1+Theta2/2),...
                RD^(1+Theta2/2),(2/(Theta2+2))*sqrt(sigma*z/CDE2S))));

            fz2=(Alpha3*(RD^2)+Alpha4*RD+c)*E-(Alpha1*RD+Alpha2)*P+RD*F;
            fm2=(Alpha3*(RD^2)+Alpha4*RD+c)*G-(Alpha1*RD+Alpha2)*Q+RD*H;
            fei2=fz2/fm2;
            
            
            fz1=Lambda*r1^(df2-df1-(Theta2-Theta1))*A-fei2*B;
            fm1=Lambda*r1^(df2-df1-(Theta2-Theta1))*C-fei2*D;
            fei1=fz1/fm1;
            
            
            ft(i,k) =  ft(i,k) + v(j) *(-fz1*beta/((z^2+z*beta)*fz1-z*fm1));
        end
        ft(i,k)=ft(i,k)*u;
        y(i,k)=-log(ft(i,k)+1)/beta;
    end
end


for k=1:1
    D_y(1,k) = tD(1) * (-3 * y(1,k) + 4 * y(2,k) - y(3,k)) / (4 * tD(2) - 3 * tD(1) - tD(3));
    D_y(2,k) = tD(2) * (y(3,k) - y(1,k)) / (tD(3) - tD(1));
    for i=3:L
        D_y(i,k) = tD(i) * (y(i - 2,k) - 4 * y(i - 1,k) + 3 * y(i,k)) / (tD(i - 2) - 4 * tD(i - 1) + 3 * tD(i));
    end
end

loglog(tD,y1(:,1),'Color',[0.75, 0, 0.75], 'Linewidth',1.5)
hold on

loglog(tD,y2(:,1),'Color',[0, 0.5, 0], 'Linewidth',1.5)
hold on

loglog(tD,y3(:,1),'Color',[0, 0, 1], 'Linewidth',1.5)
hold on

loglog(tD,y(:,1),'Color',[0.85, 0.33, 0.2], 'Linewidth',1.5)
hold on

loglog(tD,D_y1(:,1),'--','Color',[0.75, 0, 0.75], 'Linewidth',1.5)
hold on

loglog(tD,D_y2(:,1),'--','Color',[0, 0.5, 0], 'Linewidth',1.5)
hold on

loglog(tD,D_y3(:,1),'--','Color',[0, 0, 1], 'Linewidth',1.5)
hold on

loglog(tD,D_y(:,1),'--','Color',[0.85, 0.33, 0.2], 'Linewidth',1.5)
hold on

axis([10^(-2), 10^15, 10^(-2.5), 10^4.5]);
hold on

hold on
xlabel('\itt_{D}','FontName','Times New Roman','FontSize',10) 
ylabel('\itP_{wD} , P\prime_{wD}' ,'FontName','Times New Roman','FontSize',10) 

h=legend('$$\varepsilon(R_{D},t_{D})=0$$','$$\varepsilon(R_{D},t_{D})=50$$','$$\varepsilon(R_{D},t_{D})=inf$$',...
    '$$\varepsilon( {{R_D},{t_D}}) =( {{\alpha _1}{R_D} + {\alpha _2}}){t_D}+... {\alpha _3}R_D^2 + {\alpha _4}{R_D} + c$$');
set(h,'Interpreter','latex')

box off
grid on



