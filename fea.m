%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Basis truss program                              %
% Written By: Ho Yuen Henry Suen, 2021
% Content Originally from Notes and Exercise for the Course Finite Element 
% Methods 41525, Technical University of Denmark
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;close all;clc;

%--- Input file ----------------------------------------------------------%
example1                % Input file for FEM

neqn = size(X,1)*size(X,2);         % Number of equations
ne = size(IX,1);                    % Number of elements
disp(['Number of DOF ' sprintf('%d',neqn) ...
    ' Number of elements ' sprintf('%d',ne)]);

%--- Initialize arrays ---------------------------------------------------%
Kmatr=sparse(neqn,neqn);                % Stiffness matrix
P=zeros(neqn,1);                        % Force vector
D=zeros(neqn,1);                        % Displacement vector
R=zeros(neqn,1);                        % Residual vector
strain=zeros(ne,1);                     % Element strain vector
stress=zeros(ne,1);                     % Element stress vector

%--- Calculate displacements ---------------------------------------------%
[P]=buildload(X,IX,ne,P,loads,mprop);       % Build global load vector

[Kmatr]=buildstiff(X,IX,ne,mprop,Kmatr);    % Build global stiffness matrix

[Kmatr,P]=enforce(Kmatr,P,bound);           % Enforce boundary conditions

D=Kmatr\P;                                      % Solve system of equations
% D(:,1)=0;
[strain,stress]=recover(mprop,X,IX,D,ne,strain,stress); % Calculate element 
                                                        % stress and strain
                                                   
%--- Plot results --------------------------------------------------------%                                                        
PlotStructure(X,IX,ne,neqn,bound,loads,D,stress)        % Plot structure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build global load vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P]=buildload(X,IX,ne,P,loads,mprop);

for i=1:size(loads,1)
    P((loads(i,1)-1)*2+loads(i,2),1)=loads(i,3);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build global stiffness matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K]=buildstiff(X,IX,ne,mprop,K);

% This subroutine builds the global stiffness matrix from
% the local element stiffness matrices       
for e=1:ne 
    Ee=mprop(1,1);
    Ae=mprop(1,2);
    delta_x=X(IX(e,2),1)-X(IX(e,1),1);
    delta_y=X(IX(e,2),2)-X(IX(e,1),2);
    L0e=sqrt(delta_x^2+delta_y^2);
    B0=1/L0e^2*[-delta_x -delta_y delta_x delta_y]';
    k=Ee*Ae*L0e*B0*B0';
    k_number=[(IX(e,1)-1)*2+1 (IX(e,1)-1)*2+2 (IX(e,2)-1)*2+1 (IX(e,2)-1)*2+2 ];
    for i=1:4
        for j=1:4
            K(k_number(i),k_number(j))=K(k_number(i),k_number(j))+k(i,j);
        end
    end
    
end
% full(K)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Enforce boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,P]=enforce(K,P,bound);

% This subroutine enforces the support boundary conditions

for i=1:size(bound,1)
    
    K((bound(i,1)-1)*2+bound(i,2),:)=0;
    K(:,(bound(i,1)-1)*2+bound(i,2))=0;
    K((bound(i,1)-1)*2+bound(i,2),(bound(i,1)-1)*2+bound(i,2))=1;
    P((bound(i,1)-1)*2+bound(i,2))=bound(i,3);
       
end
% full(K)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calculate element strain and stress %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [strain,stress]=recover(mprop,X,IX,D,ne,strain,stress);

% This subroutine recovers the element stress, element strain, 
% and nodal reaction forces
        
for e=1:ne
    delta_x=X(IX(e,2),1)-X(IX(e,1),1);
    delta_y=X(IX(e,2),2)-X(IX(e,1),2);
    L0e=sqrt(delta_x^2+delta_y^2);
    B0=1/L0e^2*[-delta_x -delta_y delta_x delta_y];
    
    ui=D((IX(e,1)-1)*2+1)-X(IX(e,1),1);
    vi=D((IX(e,1)-1)*2+2)-X(IX(e,1),2);
    uj=D((IX(e,2)-1)*2+1)-X(IX(e,1),1);
    vj=D((IX(e,2)-1)*2+2)-X(IX(e,1),2);
    de=[ui vi uj vj]';
    
    strain(e)=B0*de;
    stress(e)=strain(e)*mprop(1,1);

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotStructure(X,IX,ne,neqn,bound,loads,D,stress)

% This subroutine plots the undeformed and deformed structure

h1=0;h2=0;
% Plotting Un-Deformed and Deformed Structure
clf
hold on
box on
for e = 1:ne
    xx = X(IX(e,1:2),1);
    yy = X(IX(e,1:2),2);
    h1=plot(xx,yy,'k:','LineWidth',1.);
    edof = [2*IX(e,1)-1 2*IX(e,1) 2*IX(e,2)-1 2*IX(e,2)];
    xx = xx + D(edof(1:2:4));
    yy = yy + D(edof(2:2:4));
    
    if stress(e)>max(stress)*10e-5&&stress(e)>0
    h2=plot(xx,yy,'b','LineWidth',3.5); 
    elseif stress(e)<-max(stress)*10e-5
    h3=plot(xx,yy,'r','LineWidth',3.5); 
    else
    h4=plot(xx,yy,'g','LineWidth',3.5);      
    end
end

plotsupports
plotloads

legend([h1 h2 h3 ],{'Undeformed state',...
                'Tension state','Compression state'})

axis equal;
hold off

end
