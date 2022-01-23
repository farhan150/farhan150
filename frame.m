%Given data%
E=20.68*(10^10)
I=2705.5*(10^-8)
A=43.87*(10^-4)
L1=2.4
L2=3.6

%Local stiffness matrix of element 1%
k_1= [A*E/L1 0 0 -A*E/L1 0 0;0 12*E*I/L1^3 6*E*I/L1^2 0 -12*E*I/L1^3 6*E*I/L1^2;0 6*E*I/L1^2 4*E*I/L1 0 -6*E*I/L1^2 2*E*I/L1;-A*E/L1 0 0 A*E/L1 0 0;0 -12*E*I/L1^3 -6*E*I/L1^2 0 12*E*I/L1^3 -6*E*I/L1^2;0 6*E*I/L1^2 2*E*I/L1 0 -6*E*I/L1^2 4*E*I/L1]

%Local or gobal stiffness matrix of element 2%
k2= [A*E/L2 0 0 -A*E/L2 0 0;0 12*E*I/L2^3 6*E*I/L2^2 0 -12*E*I/L2^3 6*E*I/L2^2;0 6*E*I/L2^2 4*E*I/L2 0 -6*E*I/L2^2 2*E*I/L2;-A*E/L2 0 0 A*E/L2 0 0;0 -12*E*I/L2^3 -6*E*I/L2^2 0 12*E*I/L2^3 -6*E*I/L2^2;0 6*E*I/L2^2 2*E*I/L2 0 -6*E*I/L2^2 4*E*I/L2]

%Local stiffness matrix of element 3%
k_3= [A*E/L1 0 0 -A*E/L1 0 0;0 12*E*I/L1^3 -6*E*I/L1^2 0 -12*E*I/L1^3 -6*E*I/L1^2;0 -6*E*I/L1^2 4*E*I/L1 0 6*E*I/L1^2 2*E*I/L1;-A*E/L1 0 0 A*E/L1 0 0;0 -12*E*I/L1^3 6*E*I/L1^2 0 12*E*I/L1^3 6*E*I/L1^2;0 -6*E*I/L1^2 2*E*I/L1 0 6*E*I/L1^2 4*E*I/L1]

%Transformation matrix of element 1 and element 3%
T=[0 1 0 0 0 0;-1 0 0 0 0 0;0 0 1 0 0 0;0 0 0 0 1 0;0 0 0 -1 0 0;0 0 0 0 0 1]

%Global stiffness matrix of element 1 and 3%
k1=(T')*(k_1)*(T)
k3=(T')*(k_3)*(T)

%degree of freedom =12 for whole structure%
dof =12

%Defining a 12 X 12 matrix assembled global stiffness matrix%
K=zeros(dof)

%Assembly of stiffness matrix in global coordinate system%
for i=1:dof
    for j=1:dof
        if i<7 && j<7
            K(i,j)=k1(i,j)
        else
            K(i,j)=K(i,j)
        end
    end
end


for i=4:9
    for j=4:9
        K(i,j)=K(i,j)+k2(i-3,j-3)
    end
end

for i=7:12
    for j=7:12
        K(i,j)=K(i,j)+k3(i-6,j-6)
    end
end

%Defining a column force matrix and filling it with known nodal forces% 
F=[0;0;0;13340;-13320;-7992;0;-13320;7992;0;0;0]

%Defining a column displacement matrix% 
d=zeros(dof,1)

%calculating nodal displcement%
d(4:9)=inv(K(4:9,4:9))*F(4:9)

%Calculating nodal forces/moment%
R=K*d

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           