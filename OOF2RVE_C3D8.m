% PROGRAM GENERATES 2D-RVE phases created by OOF2
% The RVE is created using the 3D element in ABAQUS (C3D8) when the number
% of element in thickness direction is one.
% Written by: Mehdi Shahzamanian, Mechanical Engineering Department 
% McMaster University, shahzamm@mcmaster.ca and mmshahzamanian@gmail.com
%-----------------------------------------------------------------------
clear all;
clc;
format short g
%----------------------------------------------------------------------
% SCALING PARAMETERS
%----------------------------------------------------------------------
% insert the length of the RVE in x direction 
ax=1;
% insert the length of the RVE in y direction 
ay=1;
%----------------------------------------------------------------------
az=0.1;
% insert the number of the RVE in x direction 
nx=500;
% insert the number of the RVE in y direction 
ny=500;
nz=1; 
nelm=(nx)*(ny)*(nz);
load('Matrix.mat')
load('Precipitate.mat')
xdim=size(Precipitate,1);
ydim=size(Precipitate,2);
for i=1:1:nelm;
    mat(i,1)=1;
    index=1;
    for j=1:1:ydim;
        for k=1:1:xdim;
            if (Precipitate(k,j)==i);
               mat(i,1)=2;
               idex=2;
            end         
            if (index==2)
                break
            end
        end  
        if (index==2)
            break
        end
    end
end
tic
% -----------------------------------------
% 3D NODE COORDINATES 
%------------------------------------------
display('NODES....')
nn=0;
      for k= 1:nz+1
         for j=1:ny+1
            for i=1:nx+1
               nn=nn+1;
               x(1,nn)=(i-1)*ax;
               x(2,nn)=(j-1)*ay;
               x(3,nn)=(k-1)*az;
               NODE(nn,:)=[nn,x(1,nn),x(2,nn),x(3,nn)];
            end 
         end 
      end 
      display('NODAL COORDINATES GENERATED')
%--------------------------------------------------------
%ELEMENT CONNECTIVITY TABLE FOR 8-NODE BRICK ELEMENT
%--------------------------------------------------------
display('ELEMENTS....')
n=0;
      for kk=1:(nz)
         for jj=1:(ny)
            for ii=1:(nx) 
               n=n+1;
               ELNODE(1,n,:)= ii+(jj-1)*(nx+1)+(kk-1)*(nx+1)*(ny+1);
               ELNODE(2,n,:)= (ii+1)+(jj-1)*(nx+1)+(kk-1)*(nx+1)*(ny+1);
               ELNODE(3,n,:)= (ii+1)+(jj)*(nx+1)+(kk-1)*(nx+1)*(ny+1);
               ELNODE(4,n,:)= ii+(jj)*(nx+1)+(kk-1)*(nx+1)*(ny+1);
               ELNODE(5,n,:)= ii+(jj-1)*(nx+1)+(kk)*(nx+1)*(ny+1);
               ELNODE(6,n,:)= (ii+1)+(jj-1)*(nx+1)+(kk)*(nx+1)*(ny+1);
               ELNODE(7,n,:)= (ii+1)+(jj)*(nx+1)+(kk)*(nx+1)*(ny+1);
               ELNODE(8,n,:)= ii+(jj)*(nx+1)+(kk)*(nx+1)*(ny+1);
               ELEMENTDATA(n,:)=[n, ELNODE(:,n)',mat(n,:)]; % Assign material to each element
            end
         end
      end
      nnelm=nelm;
      ncelm=nelm;
      nnc=nn; 
%----------------------------------------------------------
% To separate the elemnts at the left edge of an RVE where a
% cohesive element is to be inserted
%----------------------------------------------------------
    for i=1:nx:nelm-nx;
          if (ELEMENTDATA(i,6)~=ELEMENTDATA(i+nx,6))                       
                      jj1=ELEMENTDATA(i,5);
                      NODE(nnc+1,:)=[nnc+1,x(1,jj1),x(2,jj1),x(3,jj1)];
                      x(1,nnc+1)=x(1,jj1);
                      x(2,nnc+1)=x(2,jj1);
                      x(3,nnc+1)=x(3,jj1);
                      jj2=ELEMENTDATA(i,9);  
                      NODE(nnc+2,:)=[nnc+2,x(1,jj2),x(2,jj2),x(3,jj2)];
                      x(1,nnc+2)=x(1,jj2);
                      x(2,nnc+2)=x(2,jj2);
                      x(3,nnc+2)=x(3,jj2);
                      ELEMENTDATA(i,5)=nnc+1;
                      ELEMENTDATA(i,9)=nnc+2;
                      nnc=nnc+2;
          end
      end
%----------------------------------------------------------
% To separate the elemnts at the right edge of an RVE where a
% cohesive element is to be inserted
%----------------------------------------------------------     
     for i=nx:nx:nelm-nx;
          if (ELEMENTDATA(i,6)~=ELEMENTDATA(i+nx,6))                 
                      jj1=ELEMENTDATA(i,4);
                      NODE(nnc+1,:)=[nnc+1,x(1,jj1),x(2,jj1),x(3,jj1)];
                      x(1,nnc+1)=x(1,jj1);
                      x(2,nnc+1)=x(2,jj1);
                      x(3,nnc+1)=x(3,jj1);
                      jj2=ELEMENTDATA(i,8);  
                      NODE(nnc+2,:)=[nnc+2,x(1,jj2),x(2,jj2),x(3,jj2)];
                      x(1,nnc+2)=x(1,jj2);
                      x(2,nnc+2)=x(2,jj2);
                      x(3,nnc+2)=x(3,jj2);
                      ELEMENTDATA(i,4)=nnc+1;
                      ELEMENTDATA(i,8)=nnc+2;
                      nnc=nnc+2;
          end
      end
%----------------------------------------------------------
% To apply cohesive elements at the sides of precipitates
%----------------------------------------------------------
 for i=1:nx-1;
     if (ELEMENTDATA(i,10)~=ELEMENTDATA(i+1,10))
         if (ELEMENTDATA(i+1,10)==2)
                      jj1=ELEMENTDATA(i+1,2);  
                      NODE(nnc+1,:)=[nnc+1,x(1,jj1),x(2,jj1),x(3,jj1)];
                      x(1,nnc+1)=x(1,jj1);
                      x(2,nnc+1)=x(2,jj1);
                      x(3,nnc+1)=x(3,jj1);
                      jj2=ELEMENTDATA(i+1,5);
                      NODE(nnc+2,:)=[nnc+2,x(1,jj2),x(2,jj2),x(3,jj2)];
                      x(1,nnc+2)=x(1,jj2);
                      x(2,nnc+2)=x(2,jj2);
                      x(3,nnc+2)=x(3,jj2);
                      jj3=ELEMENTDATA(i+1,6);  
                      NODE(nnc+3,:)=[nnc+3,x(1,jj3),x(2,jj3),x(3,jj3)];
                      x(1,nnc+3)=x(1,jj3);
                      x(2,nnc+3)=x(2,jj3);
                      x(3,nnc+3)=x(3,jj3);
                      jj4=ELEMENTDATA(i+1,9);  
                      NODE(nnc+4,:)=[nnc+4,x(1,jj4),x(2,jj4),x(3,jj4)];
                      x(1,nnc+4)=x(1,jj4);
                      x(2,nnc+4)=x(2,jj4);
                      x(3,nnc+4)=x(3,jj4);
                      ELEMENTDATA(i+1,2)=nnc+1;
                      ELEMENTDATA(i+1,5)=nnc+2;
                      ELEMENTDATA(i+1,6)=nnc+3;
                      ELEMENTDATA(i+1,9)=nnc+4;
                        
                      for iii=1:nnelm;
                          if (ELEMENTDATA(iii,10)==2)
                             for jjj=2:1:9;
                                if (ELEMENTDATA(iii,jjj)==jj1)
                                    ELEMENTDATA(iii,jjj)=nnc+1;
                                elseif (ELEMENTDATA(iii,jjj)==jj2)
                                    ELEMENTDATA(iii,jjj)=nnc+2;
                                elseif (ELEMENTDATA(iii,jjj)==jj3)
                                    ELEMENTDATA(iii,jjj)=nnc+3;
                                elseif (ELEMENTDATA(iii,jjj)==jj4)
                                    ELEMENTDATA(iii,jjj)=nnc+4;
                                end
                             end
                          end
                      end
                         
                      nnc=nnc+4;
                      ELEMENTDATA(ncelm+1,1)=ncelm+1;
                      ELEMENTDATA(ncelm+1,2)=ELEMENTDATA(i+1,9);
                      ELEMENTDATA(ncelm+1,3)=ELEMENTDATA(i+1,6);
                      ELEMENTDATA(ncelm+1,4)=ELEMENTDATA(i+1,2);
                      ELEMENTDATA(ncelm+1,5)=ELEMENTDATA(i+1,5);
                      ELEMENTDATA(ncelm+1,6)=ELEMENTDATA(i,8);
                      ELEMENTDATA(ncelm+1,7)=ELEMENTDATA(i,7);
                      ELEMENTDATA(ncelm+1,8)=ELEMENTDATA(i,3);
                      ELEMENTDATA(ncelm+1,9)=ELEMENTDATA(i,4);
                      ELEMENTDATA(ncelm+1,10)=3;
                      ncelm=ncelm+1;
         elseif (ELEMENTDATA(i+1,10)==1)
                      jj1=ELEMENTDATA(i,3);  
                      NODE(nnc+1,:)=[nnc+1,x(1,jj1),x(2,jj1),x(3,jj1)];
                      x(1,nnc+1)=x(1,jj1);
                      x(2,nnc+1)=x(2,jj1);
                      x(3,nnc+1)=x(3,jj1);
                      jj2=ELEMENTDATA(i,4);
                      NODE(nnc+2,:)=[nnc+2,x(1,jj2),x(2,jj2),x(3,jj2)];
                      x(1,nnc+2)=x(1,jj2);
                      x(2,nnc+2)=x(2,jj2);
                      x(3,nnc+2)=x(3,jj2);
                      jj3=ELEMENTDATA(i,8);  
                      NODE(nnc+3,:)=[nnc+3,x(1,jj3),x(2,jj3),x(3,jj3)];
                      x(1,nnc+3)=x(1,jj3);
                      x(2,nnc+3)=x(2,jj3);
                      x(3,nnc+3)=x(3,jj3);
                      jj4=ELEMENTDATA(i,7);  
                      NODE(nnc+4,:)=[nnc+4,x(1,jj4),x(2,jj4),x(3,jj4)];
                      x(1,nnc+4)=x(1,jj4);
                      x(2,nnc+4)=x(2,jj4);
                      x(3,nnc+4)=x(3,jj4);
                      ELEMENTDATA(i,3)=nnc+1;
                      ELEMENTDATA(i,4)=nnc+2;
                      ELEMENTDATA(i,8)=nnc+3;
                      ELEMENTDATA(i,7)=nnc+4;
                      for iii=1:nnelm;
                          if (ELEMENTDATA(iii,10)==2)
                             for jjj=2:1:9;
                                if (ELEMENTDATA(iii,jjj)==jj1)
                                    ELEMENTDATA(iii,jjj)=nnc+1;
                                elseif (ELEMENTDATA(iii,jjj)==jj2)
                                    ELEMENTDATA(iii,jjj)=nnc+2;
                                elseif (ELEMENTDATA(iii,jjj)==jj3)
                                    ELEMENTDATA(iii,jjj)=nnc+3;
                                elseif (ELEMENTDATA(iii,jjj)==jj4)
                                    ELEMENTDATA(iii,jjj)=nnc+4;
                                end
                             end
                          end
                      end
                      nnc=nnc+4;
                      ELEMENTDATA(ncelm+1,1)=ncelm+1;
                      ELEMENTDATA(ncelm+1,2)=ELEMENTDATA(i+1,9);
                      ELEMENTDATA(ncelm+1,3)=ELEMENTDATA(i+1,6);
                      ELEMENTDATA(ncelm+1,4)=ELEMENTDATA(i+1,2);
                      ELEMENTDATA(ncelm+1,5)=ELEMENTDATA(i+1,5);
                      ELEMENTDATA(ncelm+1,6)=ELEMENTDATA(i,8);
                      ELEMENTDATA(ncelm+1,7)=ELEMENTDATA(i,7);
                      ELEMENTDATA(ncelm+1,8)=ELEMENTDATA(i,3);
                      ELEMENTDATA(ncelm+1,9)=ELEMENTDATA(i,4);
                      ELEMENTDATA(ncelm+1,10)=3;
                      ncelm=ncelm+1;
         end
     end
 end
 for b=2:ny;
    for i=(b-1)*nx+1:b*nx-1;

          if (ELEMENTDATA(i,10)~=ELEMENTDATA(i+1,10))
                 if (ELEMENTDATA(i+1,10)==2 && ELEMENTDATA(i+1-nx,10)==1)
                      jj1=ELEMENTDATA(i+1,2);  
                      NODE(nnc+1,:)=[nnc+1,x(1,jj1),x(2,jj1),x(3,jj1)];
                      x(1,nnc+1)=x(1,jj1);
                      x(2,nnc+1)=x(2,jj1);
                      x(3,nnc+1)=x(3,jj1);
                      jj2=ELEMENTDATA(i+1,5);
                      NODE(nnc+2,:)=[nnc+2,x(1,jj2),x(2,jj2),x(3,jj2)];
                      x(1,nnc+2)=x(1,jj2);
                      x(2,nnc+2)=x(2,jj2);
                      x(3,nnc+2)=x(3,jj2);
                      jj3=ELEMENTDATA(i+1,6);  
                      NODE(nnc+3,:)=[nnc+3,x(1,jj3),x(2,jj3),x(3,jj3)];
                      x(1,nnc+3)=x(1,jj3);
                      x(2,nnc+3)=x(2,jj3);
                      x(3,nnc+3)=x(3,jj3);
                      jj4=ELEMENTDATA(i+1,9);  
                      NODE(nnc+4,:)=[nnc+4,x(1,jj4),x(2,jj4),x(3,jj4)];
                      x(1,nnc+4)=x(1,jj4);
                      x(2,nnc+4)=x(2,jj4);
                      x(3,nnc+4)=x(3,jj4);
                      ELEMENTDATA(i+1,2)=nnc+1;
                      ELEMENTDATA(i+1,5)=nnc+2;
                      ELEMENTDATA(i+1,6)=nnc+3;
                      ELEMENTDATA(i+1,9)=nnc+4;
                      for iii=1:nnelm;
                          if (ELEMENTDATA(iii,10)==2)
                             for jjj=2:1:9;
                                if (ELEMENTDATA(iii,jjj)==jj1)
                                    ELEMENTDATA(iii,jjj)=nnc+1;
                                elseif (ELEMENTDATA(iii,jjj)==jj2)
                                    ELEMENTDATA(iii,jjj)=nnc+2;
                                elseif (ELEMENTDATA(iii,jjj)==jj3)
                                    ELEMENTDATA(iii,jjj)=nnc+3;
                                elseif (ELEMENTDATA(iii,jjj)==jj4)
                                    ELEMENTDATA(iii,jjj)=nnc+4;
                                end
                             end
                          end
                      end
                      nnc=nnc+4;
                      ELEMENTDATA(ncelm+1,1)=ncelm+1;
                      ELEMENTDATA(ncelm+1,2)=ELEMENTDATA(i+1,9);
                      ELEMENTDATA(ncelm+1,3)=ELEMENTDATA(i+1,6);
                      ELEMENTDATA(ncelm+1,4)=ELEMENTDATA(i+1,2);
                      ELEMENTDATA(ncelm+1,5)=ELEMENTDATA(i+1,5);
                      ELEMENTDATA(ncelm+1,6)=ELEMENTDATA(i,8);
                      ELEMENTDATA(ncelm+1,7)=ELEMENTDATA(i,7);
                      ELEMENTDATA(ncelm+1,8)=ELEMENTDATA(i,3);
                      ELEMENTDATA(ncelm+1,9)=ELEMENTDATA(i,4);
                      ELEMENTDATA(ncelm+1,10)=3;
                      ncelm=ncelm+1;
                  elseif (ELEMENTDATA(i+1,10)==1 && ELEMENTDATA(i-nx,10)==1)
                      jj1=ELEMENTDATA(i,3);  
                      NODE(nnc+1,:)=[nnc+1,x(1,jj1),x(2,jj1),x(3,jj1)];
                      x(1,nnc+1)=x(1,jj1);
                      x(2,nnc+1)=x(2,jj1);
                      x(3,nnc+1)=x(3,jj1);
                      jj2=ELEMENTDATA(i,4);
                      NODE(nnc+2,:)=[nnc+2,x(1,jj2),x(2,jj2),x(3,jj2)];
                      x(1,nnc+2)=x(1,jj2);
                      x(2,nnc+2)=x(2,jj2);
                      x(3,nnc+2)=x(3,jj2);
                      jj3=ELEMENTDATA(i,8);  
                      NODE(nnc+3,:)=[nnc+3,x(1,jj3),x(2,jj3),x(3,jj3)];
                      x(1,nnc+3)=x(1,jj3);
                      x(2,nnc+3)=x(2,jj3);
                      x(3,nnc+3)=x(3,jj3);
                      jj4=ELEMENTDATA(i,7);  
                      NODE(nnc+4,:)=[nnc+4,x(1,jj4),x(2,jj4),x(3,jj4)];
                      x(1,nnc+4)=x(1,jj4);
                      x(2,nnc+4)=x(2,jj4);
                      x(3,nnc+4)=x(3,jj4);
                      ELEMENTDATA(i,3)=nnc+1;
                      ELEMENTDATA(i,4)=nnc+2;
                      ELEMENTDATA(i,8)=nnc+3;
                      ELEMENTDATA(i,7)=nnc+4;
                      for iii=1:nnelm;
                          if (ELEMENTDATA(iii,10)==2)
                             for jjj=2:1:9;
                                if (ELEMENTDATA(iii,jjj)==jj1)
                                    ELEMENTDATA(iii,jjj)=nnc+1;
                                elseif (ELEMENTDATA(iii,jjj)==jj2)
                                    ELEMENTDATA(iii,jjj)=nnc+2;
                                elseif (ELEMENTDATA(iii,jjj)==jj3)
                                    ELEMENTDATA(iii,jjj)=nnc+3;
                                elseif (ELEMENTDATA(iii,jjj)==jj4)
                                    ELEMENTDATA(iii,jjj)=nnc+4;
                                end
                             end
                          end
                      end
                      nnc=nnc+4;
                      ELEMENTDATA(ncelm+1,1)=ncelm+1;
                      ELEMENTDATA(ncelm+1,2)=ELEMENTDATA(i+1,9);
                      ELEMENTDATA(ncelm+1,3)=ELEMENTDATA(i+1,6);
                      ELEMENTDATA(ncelm+1,4)=ELEMENTDATA(i+1,2);
                      ELEMENTDATA(ncelm+1,5)=ELEMENTDATA(i+1,5);
                      ELEMENTDATA(ncelm+1,6)=ELEMENTDATA(i,8);
                      ELEMENTDATA(ncelm+1,7)=ELEMENTDATA(i,7);
                      ELEMENTDATA(ncelm+1,8)=ELEMENTDATA(i,3);
                      ELEMENTDATA(ncelm+1,9)=ELEMENTDATA(i,4);
                      ELEMENTDATA(ncelm+1,10)=3;
                      ncelm=ncelm+1;
                 elseif (ELEMENTDATA(i+1,10)==2 && ELEMENTDATA(i+1-nx,10)==2 && ELEMENTDATA(i-nx,10)==1)
                      jj1=ELEMENTDATA(i+1,5);
                      NODE(nnc+1,:)=[nnc+1,x(1,jj1),x(2,jj1),x(3,jj1)];
                      x(1,nnc+1)=x(1,jj1);
                      x(2,nnc+1)=x(2,jj1);
                      x(3,nnc+1)=x(3,jj1);
                      jj2=ELEMENTDATA(i+1,9);  
                      NODE(nnc+2,:)=[nnc+2,x(1,jj2),x(2,jj2),x(3,jj2)];
                      x(1,nnc+2)=x(1,jj2);
                      x(2,nnc+2)=x(2,jj2);
                      x(3,nnc+2)=x(3,jj2);
                      ELEMENTDATA(i+1,5)=nnc+1;
                      ELEMENTDATA(i+1,9)=nnc+2;
                      for iii=1:nnelm;
                          if (ELEMENTDATA(iii,10)==2)
                             for jjj=2:1:9;
                                if (ELEMENTDATA(iii,jjj)==jj1)
                                    ELEMENTDATA(iii,jjj)=nnc+1;
                                elseif (ELEMENTDATA(iii,jjj)==jj2)
                                    ELEMENTDATA(iii,jjj)=nnc+2;
                                end
                             end
                          end
                      end
                      nnc=nnc+2;
                      ELEMENTDATA(ncelm+1,1)=ncelm+1;
                      ELEMENTDATA(ncelm+1,2)=ELEMENTDATA(i+1,9);
                      ELEMENTDATA(ncelm+1,3)=ELEMENTDATA(i+1,6);
                      ELEMENTDATA(ncelm+1,4)=ELEMENTDATA(i+1,2);
                      ELEMENTDATA(ncelm+1,5)=ELEMENTDATA(i+1,5);
                      ELEMENTDATA(ncelm+1,6)=ELEMENTDATA(i,8);
                      ELEMENTDATA(ncelm+1,7)=ELEMENTDATA(i,7);
                      ELEMENTDATA(ncelm+1,8)=ELEMENTDATA(i,3);
                      ELEMENTDATA(ncelm+1,9)=ELEMENTDATA(i,4);
                      ELEMENTDATA(ncelm+1,10)=3;
                      ncelm=ncelm+1;
                elseif (ELEMENTDATA(i+1,10)==2 && ELEMENTDATA(i+1-nx,10)==2 && ELEMENTDATA(i-nx,10)==2)
                      jj1=ELEMENTDATA(i+1,2);  
                      NODE(nnc+1,:)=[nnc+1,x(1,jj1),x(2,jj1),x(3,jj1)];
                      x(1,nnc+1)=x(1,jj1);
                      x(2,nnc+1)=x(2,jj1);
                      x(3,nnc+1)=x(3,jj1);
                      jj2=ELEMENTDATA(i+1,5);
                      NODE(nnc+2,:)=[nnc+2,x(1,jj2),x(2,jj2),x(3,jj2)];
                      x(1,nnc+2)=x(1,jj2);
                      x(2,nnc+2)=x(2,jj2);
                      x(3,nnc+2)=x(3,jj2);
                      jj3=ELEMENTDATA(i+1,6);  
                      NODE(nnc+3,:)=[nnc+3,x(1,jj3),x(2,jj3),x(3,jj3)];
                      x(1,nnc+3)=x(1,jj3);
                      x(2,nnc+3)=x(2,jj3);
                      x(3,nnc+3)=x(3,jj3);
                      jj4=ELEMENTDATA(i+1,9);  
                      NODE(nnc+4,:)=[nnc+4,x(1,jj4),x(2,jj4),x(3,jj4)];
                      x(1,nnc+4)=x(1,jj4);
                      x(2,nnc+4)=x(2,jj4);
                      x(3,nnc+4)=x(3,jj4);
                      ELEMENTDATA(i+1,2)=nnc+1;
                      ELEMENTDATA(i+1,5)=nnc+2;
                      ELEMENTDATA(i+1,6)=nnc+3;
                      ELEMENTDATA(i+1,9)=nnc+4;
                      for iii=1:nnelm;
                          if (ELEMENTDATA(iii,10)==2)
                             for jjj=2:1:9;
                                if (ELEMENTDATA(iii,jjj)==jj1)
                                    ELEMENTDATA(iii,jjj)=nnc+1;
                                elseif (ELEMENTDATA(iii,jjj)==jj2)
                                    ELEMENTDATA(iii,jjj)=nnc+2;
                                elseif (ELEMENTDATA(iii,jjj)==jj3)
                                    ELEMENTDATA(iii,jjj)=nnc+3;
                                elseif (ELEMENTDATA(iii,jjj)==jj4)
                                    ELEMENTDATA(iii,jjj)=nnc+4;
                                end
                             end
                          end
                      end
                      nnc=nnc+4;
                      ELEMENTDATA(ncelm+1,1)=ncelm+1;
                      ELEMENTDATA(ncelm+1,2)=ELEMENTDATA(i+1,9);
                      ELEMENTDATA(ncelm+1,3)=ELEMENTDATA(i+1,6);
                      ELEMENTDATA(ncelm+1,4)=ELEMENTDATA(i+1,2);
                      ELEMENTDATA(ncelm+1,5)=ELEMENTDATA(i+1,5);
                      ELEMENTDATA(ncelm+1,6)=ELEMENTDATA(i,8);
                      ELEMENTDATA(ncelm+1,7)=ELEMENTDATA(i,7);
                      ELEMENTDATA(ncelm+1,8)=ELEMENTDATA(i,3);
                      ELEMENTDATA(ncelm+1,9)=ELEMENTDATA(i,4);
                      ELEMENTDATA(ncelm+1,10)=3;
                      ncelm=ncelm+1;  
                           
                elseif (ELEMENTDATA(i+1,10)==1 && ELEMENTDATA(i-nx,10)==2 &&  ELEMENTDATA(i+1-nx,10)==1)
                      jj1=ELEMENTDATA(i,4);
                      NODE(nnc+1,:)=[nnc+1,x(1,jj1),x(2,jj1),x(3,jj1)];
                      x(1,nnc+1)=x(1,jj1);
                      x(2,nnc+1)=x(2,jj1);
                      x(3,nnc+1)=x(3,jj1);
                      jj2=ELEMENTDATA(i,8);  
                      NODE(nnc+2,:)=[nnc+2,x(1,jj2),x(2,jj2),x(3,jj2)];
                      x(1,nnc+2)=x(1,jj2);
                      x(2,nnc+2)=x(2,jj2);
                      x(3,nnc+2)=x(3,jj2);
                      ELEMENTDATA(i,4)=nnc+1;
                      ELEMENTDATA(i,8)=nnc+2;
                      for iii=1:nnelm;
                          if (ELEMENTDATA(iii,10)==2)
                             for jjj=2:1:9;
                                if (ELEMENTDATA(iii,jjj)==jj1)
                                    ELEMENTDATA(iii,jjj)=nnc+1;
                                elseif (ELEMENTDATA(iii,jjj)==jj2)
                                    ELEMENTDATA(iii,jjj)=nnc+2;
                                end
                             end
                          end
                      end
                      nnc=nnc+2;
                      ELEMENTDATA(ncelm+1,1)=ncelm+1;
                      ELEMENTDATA(ncelm+1,2)=ELEMENTDATA(i+1,9);
                      ELEMENTDATA(ncelm+1,3)=ELEMENTDATA(i+1,6);
                      ELEMENTDATA(ncelm+1,4)=ELEMENTDATA(i+1,2);
                      ELEMENTDATA(ncelm+1,5)=ELEMENTDATA(i+1,5);
                      ELEMENTDATA(ncelm+1,6)=ELEMENTDATA(i,8);
                      ELEMENTDATA(ncelm+1,7)=ELEMENTDATA(i,7);
                      ELEMENTDATA(ncelm+1,8)=ELEMENTDATA(i,3);
                      ELEMENTDATA(ncelm+1,9)=ELEMENTDATA(i,4);
                      ELEMENTDATA(ncelm+1,10)=3;
                      ncelm=ncelm+1;
                
                 elseif (ELEMENTDATA(i+1,10)==1 && ELEMENTDATA(i-nx,10)==2 &&  ELEMENTDATA(i+1-nx,10)==2)
                      jj1=ELEMENTDATA(i,3);  
                      NODE(nnc+1,:)=[nnc+1,x(1,jj1),x(2,jj1),x(3,jj1)];
                      x(1,nnc+1)=x(1,jj1);
                      x(2,nnc+1)=x(2,jj1);
                      x(3,nnc+1)=x(3,jj1);
                      jj2=ELEMENTDATA(i,4);
                      NODE(nnc+2,:)=[nnc+2,x(1,jj2),x(2,jj2),x(3,jj2)];
                      x(1,nnc+2)=x(1,jj2);
                      x(2,nnc+2)=x(2,jj2);
                      x(3,nnc+2)=x(3,jj2);
                      jj3=ELEMENTDATA(i,8);  
                      NODE(nnc+3,:)=[nnc+3,x(1,jj3),x(2,jj3),x(3,jj3)];
                      x(1,nnc+3)=x(1,jj3);
                      x(2,nnc+3)=x(2,jj3);
                      x(3,nnc+3)=x(3,jj3);
                      jj4=ELEMENTDATA(i,7);  
                      NODE(nnc+4,:)=[nnc+4,x(1,jj4),x(2,jj4),x(3,jj4)];
                      x(1,nnc+4)=x(1,jj4);
                      x(2,nnc+4)=x(2,jj4);
                      x(3,nnc+4)=x(3,jj4);
                      ELEMENTDATA(i,3)=nnc+1;
                      ELEMENTDATA(i,4)=nnc+2;
                      ELEMENTDATA(i,8)=nnc+3;
                      ELEMENTDATA(i,7)=nnc+4;
                      for iii=1:nnelm;
                          if (ELEMENTDATA(iii,10)==2)
                             for jjj=2:1:9;
                                if (ELEMENTDATA(iii,jjj)==jj1)
                                    ELEMENTDATA(iii,jjj)=nnc+1;
                                elseif (ELEMENTDATA(iii,jjj)==jj2)
                                    ELEMENTDATA(iii,jjj)=nnc+2;
                                elseif (ELEMENTDATA(iii,jjj)==jj3)
                                    ELEMENTDATA(iii,jjj)=nnc+3;
                                elseif (ELEMENTDATA(iii,jjj)==jj4)
                                    ELEMENTDATA(iii,jjj)=nnc+4;
                                end
                             end
                          end
                      end
                      nnc=nnc+4;
                      ELEMENTDATA(ncelm+1,1)=ncelm+1;
                      ELEMENTDATA(ncelm+1,2)=ELEMENTDATA(i+1,9);
                      ELEMENTDATA(ncelm+1,3)=ELEMENTDATA(i+1,6);
                      ELEMENTDATA(ncelm+1,4)=ELEMENTDATA(i+1,2);
                      ELEMENTDATA(ncelm+1,5)=ELEMENTDATA(i+1,5);
                      ELEMENTDATA(ncelm+1,6)=ELEMENTDATA(i,8);
                      ELEMENTDATA(ncelm+1,7)=ELEMENTDATA(i,7);
                      ELEMENTDATA(ncelm+1,8)=ELEMENTDATA(i,3);
                      ELEMENTDATA(ncelm+1,9)=ELEMENTDATA(i,4);
                      ELEMENTDATA(ncelm+1,10)=3;
                      ncelm=ncelm+1;
                      
                 elseif (ELEMENTDATA(i+1,10)==1 && ELEMENTDATA(i-nx,10)==2)
                      jj1=ELEMENTDATA(i,4);
                      NODE(nnc+1,:)=[nnc+1,x(1,jj1),x(2,jj1),x(3,jj1)];
                      x(1,nnc+1)=x(1,jj1);
                      x(2,nnc+1)=x(2,jj1);
                      x(3,nnc+1)=x(3,jj1);
                      jj2=ELEMENTDATA(i,8);  
                      NODE(nnc+2,:)=[nnc+2,x(1,jj2),x(2,jj2),x(3,jj2)];
                      x(1,nnc+2)=x(1,jj2);
                      x(2,nnc+2)=x(2,jj2);
                      x(3,nnc+2)=x(3,jj2);
                      ELEMENTDATA(i,4)=nnc+1;
                      ELEMENTDATA(i,8)=nnc+2;
                      for iii=1:nnelm;
                          if (ELEMENTDATA(iii,10)==2)
                             for jjj=2:1:9;
                                if (ELEMENTDATA(iii,jjj)==jj1)
                                    ELEMENTDATA(iii,jjj)=nnc+1;
                                elseif (ELEMENTDATA(iii,jjj)==jj2)
                                    ELEMENTDATA(iii,jjj)=nnc+2;
                                end
                             end
                          end
                      end
                      nnc=nnc+2;
                      ELEMENTDATA(ncelm+1,1)=ncelm+1;
                      ELEMENTDATA(ncelm+1,2)=ELEMENTDATA(i+1,9);
                      ELEMENTDATA(ncelm+1,3)=ELEMENTDATA(i+1,6);
                      ELEMENTDATA(ncelm+1,4)=ELEMENTDATA(i+1,2);
                      ELEMENTDATA(ncelm+1,5)=ELEMENTDATA(i+1,5);
                      ELEMENTDATA(ncelm+1,6)=ELEMENTDATA(i,8);
                      ELEMENTDATA(ncelm+1,7)=ELEMENTDATA(i,7);
                      ELEMENTDATA(ncelm+1,8)=ELEMENTDATA(i,3);
                      ELEMENTDATA(ncelm+1,9)=ELEMENTDATA(i,4);
                      ELEMENTDATA(ncelm+1,10)=3;
                      ncelm=ncelm+1;     
                 end 
          end
    end
end    
    
%-----------------------------------------------------------------
% To apply cohesive elements at the bottom sides of precipitates
%-----------------------------------------------------------------
      for i=nx+1:nelm; 
                  if (ELEMENTDATA(i,10)==2 && ELEMENTDATA(i,10)~=ELEMENTDATA(i-nx,10) && ELEMENTDATA(i,10)==ELEMENTDATA(i+1,10) && ELEMENTDATA(i+1-nx,10)~=2 && rem(i,nx)~=0) 
                      jj1=ELEMENTDATA(i,3);
                      NODE(nnc+1,:)=[nnc+1,x(1,jj1),x(2,jj1),x(3,jj1)];
                      x(1,nnc+1)=x(1,jj1);
                      x(2,nnc+1)=x(2,jj1);
                      x(3,nnc+1)=x(3,jj1);
                      jj2=ELEMENTDATA(i,7);  
                      NODE(nnc+2,:)=[nnc+2,x(1,jj2),x(2,jj2),x(3,jj2)];
                      x(1,nnc+2)=x(1,jj2);
                      x(2,nnc+2)=x(2,jj2);
                      x(3,nnc+2)=x(3,jj2);
                      ELEMENTDATA(i,3)=nnc+1;
                      ELEMENTDATA(i,7)=nnc+2;
                      for iii=1:nnelm;
                          if (ELEMENTDATA(iii,10)==2)
                             for jjj=2:1:9;
                                if (ELEMENTDATA(iii,jjj)==jj1)
                                    ELEMENTDATA(iii,jjj)=nnc+1;
                                elseif (ELEMENTDATA(iii,jjj)==jj2)
                                    ELEMENTDATA(iii,jjj)=nnc+2;
                                end
                             end
                          end
                      end
                      nnc=nnc+2;
                      ELEMENTDATA(ncelm+1,1)=ncelm+1;
                      ELEMENTDATA(ncelm+1,2)=ELEMENTDATA(i-nx,5);
                      ELEMENTDATA(ncelm+1,3)=ELEMENTDATA(i-nx,4);
                      ELEMENTDATA(ncelm+1,4)=ELEMENTDATA(i-nx,8);
                      ELEMENTDATA(ncelm+1,5)=ELEMENTDATA(i-nx,9);
                      ELEMENTDATA(ncelm+1,6)=ELEMENTDATA(i,2);
                      ELEMENTDATA(ncelm+1,7)=ELEMENTDATA(i,3);
                      ELEMENTDATA(ncelm+1,8)=ELEMENTDATA(i,7);
                      ELEMENTDATA(ncelm+1,9)=ELEMENTDATA(i,6);
                      ELEMENTDATA(ncelm+1,10)=3;
                      ncelm=ncelm+1; 
                elseif (ELEMENTDATA(i,10)==2 && ELEMENTDATA(i,10)~=ELEMENTDATA(i-nx,10) && ELEMENTDATA(i,10)~=ELEMENTDATA(i+1,10) && rem(i,nx)~=0) 
                      ELEMENTDATA(ncelm+1,1)=ncelm+1;
                      ELEMENTDATA(ncelm+1,2)=ELEMENTDATA(i-nx,5);
                      ELEMENTDATA(ncelm+1,3)=ELEMENTDATA(i-nx,4);
                      ELEMENTDATA(ncelm+1,4)=ELEMENTDATA(i-nx,8);
                      ELEMENTDATA(ncelm+1,5)=ELEMENTDATA(i-nx,9);
                      ELEMENTDATA(ncelm+1,6)=ELEMENTDATA(i,2);
                      ELEMENTDATA(ncelm+1,7)=ELEMENTDATA(i,3);
                      ELEMENTDATA(ncelm+1,8)=ELEMENTDATA(i,7);
                      ELEMENTDATA(ncelm+1,9)=ELEMENTDATA(i,6);
                      ELEMENTDATA(ncelm+1,10)=3;
                      ncelm=ncelm+1;
                 elseif (ELEMENTDATA(i,10)==2 && ELEMENTDATA(i,10)~=ELEMENTDATA(i-nx,10) && ELEMENTDATA(i,10)==ELEMENTDATA(i+1,10) && ELEMENTDATA(i+1-nx,10)==2 && rem(i,nx)~=0)      
                      ELEMENTDATA(ncelm+1,1)=ncelm+1;
                      ELEMENTDATA(ncelm+1,2)=ELEMENTDATA(i-nx,5);
                      ELEMENTDATA(ncelm+1,3)=ELEMENTDATA(i-nx,4);
                      ELEMENTDATA(ncelm+1,4)=ELEMENTDATA(i-nx,8);
                      ELEMENTDATA(ncelm+1,5)=ELEMENTDATA(i-nx,9);
                      ELEMENTDATA(ncelm+1,6)=ELEMENTDATA(i,2);
                      ELEMENTDATA(ncelm+1,7)=ELEMENTDATA(i,3);
                      ELEMENTDATA(ncelm+1,8)=ELEMENTDATA(i,7);
                      ELEMENTDATA(ncelm+1,9)=ELEMENTDATA(i,6);
                      ELEMENTDATA(ncelm+1,10)=3; 
                      ncelm=ncelm+1;
                elseif (ELEMENTDATA(i,10)==2 && ELEMENTDATA(i,10)~=ELEMENTDATA(i-nx,10) && rem(i,nx)==0)      
                      ELEMENTDATA(ncelm+1,1)=ncelm+1;
                      ELEMENTDATA(ncelm+1,2)=ELEMENTDATA(i-nx,5);
                      ELEMENTDATA(ncelm+1,3)=ELEMENTDATA(i-nx,4);
                      ELEMENTDATA(ncelm+1,4)=ELEMENTDATA(i-nx,8);
                      ELEMENTDATA(ncelm+1,5)=ELEMENTDATA(i-nx,9);
                      ELEMENTDATA(ncelm+1,6)=ELEMENTDATA(i,2);
                      ELEMENTDATA(ncelm+1,7)=ELEMENTDATA(i,3);
                      ELEMENTDATA(ncelm+1,8)=ELEMENTDATA(i,7);
                      ELEMENTDATA(ncelm+1,9)=ELEMENTDATA(i,6);
                      ELEMENTDATA(ncelm+1,10)=3; 
                      ncelm=ncelm+1;
                end
      
      end
%-----------------------------------------------------------------
% To apply cohesive elements at the top sides of precipitates
%-----------------------------------------------------------------      
      for i=1:(ny-1)*nx;
               if (ELEMENTDATA(i,10)==2 && ELEMENTDATA(i,10)~=ELEMENTDATA(i+nx,10) && ELEMENTDATA(i,10)==ELEMENTDATA(i+1,10) && ELEMENTDATA(i+1+nx,10)==1 && rem(i,nx)~=0) 
                      jj1=ELEMENTDATA(i,4);
                      NODE(nnc+1,:)=[nnc+1,x(1,jj1),x(2,jj1),x(3,jj1)];
                      x(1,nnc+1)=x(1,jj1);
                      x(2,nnc+1)=x(2,jj1);
                      x(3,nnc+1)=x(3,jj1);
                      jj2=ELEMENTDATA(i,8);  
                      NODE(nnc+2,:)=[nnc+2,x(1,jj2),x(2,jj2),x(3,jj2)];
                      x(1,nnc+2)=x(1,jj2);
                      x(2,nnc+2)=x(2,jj2);
                      x(3,nnc+2)=x(3,jj2);
                      ELEMENTDATA(i,4)=nnc+1;
                      ELEMENTDATA(i,8)=nnc+2;
                      for iii=1:nnelm;
                          if (ELEMENTDATA(iii,10)==2)
                             for jjj=2:1:9;
                                if (ELEMENTDATA(iii,jjj)==jj1)
                                    ELEMENTDATA(iii,jjj)=nnc+1;
                                elseif (ELEMENTDATA(iii,jjj)==jj2)
                                    ELEMENTDATA(iii,jjj)=nnc+2;
                                end
                             end
                          end
                      end
                      nnc=nnc+2;
                      ELEMENTDATA(ncelm+1,1)=ncelm+1;
                      ELEMENTDATA(ncelm+1,2)=ELEMENTDATA(i+nx,2);
                      ELEMENTDATA(ncelm+1,3)=ELEMENTDATA(i+nx,3);
                      ELEMENTDATA(ncelm+1,4)=ELEMENTDATA(i+nx,7);
                      ELEMENTDATA(ncelm+1,5)=ELEMENTDATA(i+nx,6);
                      ELEMENTDATA(ncelm+1,6)=ELEMENTDATA(i,5);
                      ELEMENTDATA(ncelm+1,7)=ELEMENTDATA(i,4);
                      ELEMENTDATA(ncelm+1,8)=ELEMENTDATA(i,8);
                      ELEMENTDATA(ncelm+1,9)=ELEMENTDATA(i,9);
                      ELEMENTDATA(ncelm+1,10)=3;
                      ncelm=ncelm+1; 
                elseif (ELEMENTDATA(i,10)==2 && ELEMENTDATA(i,10)~=ELEMENTDATA(i+nx,10) && ELEMENTDATA(i,10)~=ELEMENTDATA(i+1,10) && rem(i,nx)~=0) 
                      ELEMENTDATA(ncelm+1,1)=ncelm+1;
                      ELEMENTDATA(ncelm+1,2)=ELEMENTDATA(i+nx,2);
                      ELEMENTDATA(ncelm+1,3)=ELEMENTDATA(i+nx,3);
                      ELEMENTDATA(ncelm+1,4)=ELEMENTDATA(i+nx,7);
                      ELEMENTDATA(ncelm+1,5)=ELEMENTDATA(i+nx,6);
                      ELEMENTDATA(ncelm+1,6)=ELEMENTDATA(i,5);
                      ELEMENTDATA(ncelm+1,7)=ELEMENTDATA(i,4);
                      ELEMENTDATA(ncelm+1,8)=ELEMENTDATA(i,8);
                      ELEMENTDATA(ncelm+1,9)=ELEMENTDATA(i,9);
                      ELEMENTDATA(ncelm+1,10)=3;
                      ncelm=ncelm+1;
                elseif (ELEMENTDATA(i,10)==2 && ELEMENTDATA(i,10)~=ELEMENTDATA(i+nx,10) && ELEMENTDATA(i,10)==ELEMENTDATA(i+1,10) && ELEMENTDATA(i+1+nx,10)==2 && rem(i,nx)~=0) 
                      ELEMENTDATA(ncelm+1,1)=ncelm+1;
                      ELEMENTDATA(ncelm+1,2)=ELEMENTDATA(i+nx,2);
                      ELEMENTDATA(ncelm+1,3)=ELEMENTDATA(i+nx,3);
                      ELEMENTDATA(ncelm+1,4)=ELEMENTDATA(i+nx,7);
                      ELEMENTDATA(ncelm+1,5)=ELEMENTDATA(i+nx,6);
                      ELEMENTDATA(ncelm+1,6)=ELEMENTDATA(i,5);
                      ELEMENTDATA(ncelm+1,7)=ELEMENTDATA(i,4);
                      ELEMENTDATA(ncelm+1,8)=ELEMENTDATA(i,8);
                      ELEMENTDATA(ncelm+1,9)=ELEMENTDATA(i,9);
                      ELEMENTDATA(ncelm+1,10)=3;
                      ncelm=ncelm+1;
                elseif (ELEMENTDATA(i,10)==2 && ELEMENTDATA(i,10)~=ELEMENTDATA(i+nx,10) && rem(i,nx)==0) 
                      ELEMENTDATA(ncelm+1,1)=ncelm+1;
                      ELEMENTDATA(ncelm+1,2)=ELEMENTDATA(i+nx,2);
                      ELEMENTDATA(ncelm+1,3)=ELEMENTDATA(i+nx,3);
                      ELEMENTDATA(ncelm+1,4)=ELEMENTDATA(i+nx,7);
                      ELEMENTDATA(ncelm+1,5)=ELEMENTDATA(i+nx,6);
                      ELEMENTDATA(ncelm+1,6)=ELEMENTDATA(i,5);
                      ELEMENTDATA(ncelm+1,7)=ELEMENTDATA(i,4);
                      ELEMENTDATA(ncelm+1,8)=ELEMENTDATA(i,8);
                      ELEMENTDATA(ncelm+1,9)=ELEMENTDATA(i,9);
                      ELEMENTDATA(ncelm+1,10)=3;
                      ncelm=ncelm+1;
                end
      end
%--------------------------------------------------------------------------      
      nelm=ncelm;
      nn=nnc;
ELEMENTCONN=ELEMENTDATA(:,1:9);
display('ELEMENT CONNECTIVITY TABLE GENERATED')
%--------------------------------------------------------------------------
% OPEN A FILE FOR WRITING THE *INP FILE
%--------------------------------------------------------------------------
fid = fopen('OOF2RVE_C3D8.inp', 'w'); % Name of ABAQUS Input File
% PRINT FILE HEADER
fprintf(fid, '*HEADING \n');
fprintf(fid, 'RVE \n');
% PRINT NODE NUMBERS AND COORDINATES
fprintf(fid,'*NODE \n');
fprintf(fid, '% i,   % f,   % f,   % f \n', NODE')
display('NODAL COORDINATES GENERATED')
for ijk=1:nelm
    if (ijk<=nnelm)
    ELEMENTCONNC3D8(ijk,1:9)=ELEMENTCONN(ijk,1:9);
    else
    ELEMENTCONNCOH3D8(ijk-nnelm,1:9)=ELEMENTCONN(ijk,1:9);
    end
end
fprintf(fid,'*ELEMENT, TYPE=C3D8\n');
fprintf(fid, '%i,   %i,   %i,   %i,   %i,   %i,   %i,   %i,   %i,\n', ELEMENTCONNC3D8')
fprintf(fid,'*ELEMENT, TYPE=COH3D8\n');
fprintf(fid, '%i,   %i,   %i,   %i,   %i,   %i,   %i,   %i,   %i,\n', ELEMENTCONNCOH3D8')
%----------------------------------------------
% SORT ELEMENTS BASED ON MATERIALS
%------------------------------------------------
display('GENERATING ELEMENT SETS BASED ON MATERIAL TYPES')
BB = sortrows(ELEMENTDATA,10);%sort elements based on material type
MATERIALS=unique(BB(:,10));% Material identifiers
NMAT=size(MATERIALS,1);%Total number of materials
% Separate the elements into element sets based on the material type =ELEMSET
for jjj=1:NMAT
    for iii=1:nelm
        if ELEMENTDATA(iii,10)==MATERIALS(jjj)
            ELEMSET(iii,jjj)=ELEMENTDATA(iii,1);
        end 
    end
end
ELEMSET;
%-----------------------------------------------------------
%GENERATE ABAQUS INPUT FILE
%-----------------------------------------------------------
% DEFINE ELEMENT SETS BASED ON MATERIAL TYPES
for i=1:NMAT
    ELEMSETSORT=nonzeros(ELEMSET(:,i));    
    fprintf(fid,'\n*ELSET,ELSET=% c \n',i+65);
    fprintf(fid,'% i,  % i,  % i, % i, % i,  % i,  % i,  % i,  % i,  % i, % i, % i,  % i,  % i,  % i \n', ELEMSETSORT');
end
clear ELEMSET;
clear ELEMSETSORT;
%---------------------------------------------------------------------------------------------
% DEFINE SOLID SECTIONS BASED ON MATERIAL TYPES
for i=1:NMAT-1
    fprintf(fid,'\n*SOLID SECTION, ELSET=% c,  MATERIAL=% c % c \n',i+65, i+65, i+65);
end
i=NMAT;
fprintf(fid,'\n*Cohesive SECTION, ELSET=% c,  MATERIAL=% c % c, , response=TRACTION SEPARATION  \n',i+65, i+65, i+65);
% DEFINE ELASTIC MATERIAL PROPERTIES (Modulus, poissons's ratio)
MODULUS=10000*ones(NMAT,1);
POISSONRATIO=0.3*ones(NMAT,1); 
MATPROPELAST=[MODULUS POISSONRATIO];
for i=1:NMAT  
    fprintf(fid,'\n*MATERIAL, NAME=% c % c \n',i+65, i+65);
    fprintf(fid,'*ELASTIC \n');
    fprintf(fid,'% f,  % f \n',MATPROPELAST(i,1),MATPROPELAST(i,2));
end
display('DATA WRITTEN TO ABAQUS *inp FILE')
toc



