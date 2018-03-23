# Fast marching
# input: F=speed function, Seeds=2 lines, n columns matrix, 
# output: T=matrix with time that each pixel was hit by the front
#-----------------------------------
# There are 3 pixel status in the matrix Status (has dim(F)):
#   -1 = alive (processed)
#    0 = narrow band (trial = pending pixels to have their distance inspected, process that result in finding the smallest distance)
#    Z+ far (any integer >0 means not yet used and index matrix T)

msfm2d4 = function(F, Seeds)
{
require(PolynomF)
eps = 1.110223e-16;
MAX_timer = 40000;
alive = -1
far = 0
#narrowband is any positive integer, also index of T

Status = matrix(0,dim(F)[1],dim(F)[2]);
storage.mode(Status) <-"integer"

T = matrix(Inf,dim(F)[1],dim(F)[2]); #initial value is +infinite
storage.mode(Status) <-"integer" #it will change to double automatically just in case; performance issues?

# List substituing the heap - improve this part for O(nlogn)!
neigh_free = 100000;
nLtrial=0;
Ltrial = matrix(0,4,neigh_free);
#storage.mode(Ltrial) <-"integer"
timer = 0;

# Neighbours
neigh = matrix(c(-1, 0, 1, 0, 0, -1, 0, 1), 4, 2, byrow=T);
storage.mode(neigh)<-"integer" #saves a lot of memory

# Part (1) all seeds are Alive------------------------------------------<
for (z in 1:(dim(Seeds)[2])){
    # seed pixel
    x= Seeds[1,z]; y=Seeds[2,z];    
    Status[x,y]=alive; T[x,y]=0;
}#----------------------------------------------------------------------------------->

#Part (2) all neighbors of Alive are in Narrow Band (add them to list)------------------------------------------<
for (z in 1:(dim(Seeds)[2])){
    # seed pixel
    x= Seeds[1,z]; y=Seeds[2,z];
    for (k in 1:4) {
        # Location of neighbour
        i=x+neigh[k,1]; j=y+neigh[k,2];
        # Is index valid and the pixel is in far, after all you don't update whoever is already in the list 
        if( (i>0)&&(j>0)&&(i<=dim(F)[1])&&(j<=dim(F)[2])&&(Status[i,j]==far) )
	{
            Tt=(1/(F[i,j]+eps)); #temporary time calculation
            # Update distance in neigbour list or add to neigbour list
            #unnecessaryif(is.finite(Status[i,j]))#is it in the trial list?
            #    Ltrial[1,Status[i,j]] = min(Tt,Ltrial[1,Status[i,j]]) #yes, it is
            #else{
            nLtrial=nLtrial+1;# no, it isn't, so make it part of the trial list
            # If running out of memory at a new block
            if(nLtrial>neigh_free){ neigh_free = neigh_free +100000; Ltrial[1,neigh_free]=0;}
            Ltrial[,nLtrial]=c(Tt,i,j,timer);
            Status[i,j]=nLtrial;      #}
	    T[i,j]=Ltrial[1,Status[i,j]]
        }
    }
}#end initialization here------------------------------------------------------------------------------------>

#Part (3) start iterative process of propagating boundary ---------------------------------------------------<
#Loop through ALL pixels of the image
for (itt in 1:length(F)){  
    timer = timer + 1
    if (nLtrial==0) break  
    # Get the pixel from narrow list with smallest distance value and set it to current pixel location
    t = min(Ltrial[1,1:nLtrial]);
    index = which(Ltrial[1,1:nLtrial] == t, arr.ind=T);
    index = index[1]; #can return more than one value	
    x=Ltrial[2,index]; y=Ltrial[3,index];
    Status[x,y]=alive;
    T[x,y]=Ltrial[1,index]; #receive the min value
    
    # Remove min value by replacing it with the last value in the array
    if(index<nLtrial){
        Ltrial[,index]=Ltrial[,nLtrial];
	Ltrial[,nLtrial] = 0
        x2=Ltrial[2,index]; y2=Ltrial[3,index];
        Status[x2,y2]=index; 
    }
    nLtrial = nLtrial-1;
    
    if( (timer - Ltrial[4,index]) < MAX_timer ){    
    # Loop through all 4 neighbours of current pixel to update their time
	    for (k in 1:4){	
	        # Location of neighbour
	        i=x+neigh[k,1]; j=y+neigh[k,2];
	        
	        # Is index valid and the pixel not yet Alive?
	        if( (i>0)&&(j>0) && (i<=dim(F)[1]) && (j<=dim(F)[2]) && (Status[i,j]!=alive) ){	            
	            # Boundary and Alive check -> current patch
	            Tpatch=matrix(Inf,3,3); 
	            for (nx in -1:1)
	                for (ny in -1:1){
	                    ii=i+nx; jj=j+ny;			
	                    if((ii>0)&&(jj>0)&&(ii<=dim(F)[1])&&(jj<=dim(F)[2])&&(is.finite(T[ii,jj]))){
				#print(paste(ii,jj,sep=','))
	                        Tpatch[nx+2,ny+2]=T[ii,jj];                                        
				}
	            	}
	            # The values in order is 0 if no neighbours in that direction
	            # 1 if 1e order derivatives is used 
	            Order=matrix(0,1,2);
	            Tm=array(0,2)
	            # Make 1e order derivatives in x and y direction
	            Tm[1] = min( Tpatch[1,2] , Tpatch[3,2]); if(is.finite(Tm[1])) Order[1]=1; 
	            Tm[2] = min( Tpatch[2,1] , Tpatch[2,3]); if(is.finite(Tm[2])) Order[2]=1; 
	                        
	            # Calculate the distance using x and y direction
	            Coeff = c(0, 0, -1/(F[i,j]^2+eps));	    
		    for (ind in 1:2){
	              if(Order[ind])
	                  Coeff=Coeff+c(1, -2*Tm[ind], Tm[ind]^2); 
		    }
		    #if (i==49&j==53) 		browser()
		    v = polynom()
		    p = Coeff[1]*v^2 + Coeff[2]*v + Coeff[3]
		    Tt=NULL
	            Tt=solve(p) 		    
		    if (!is.complex(Tt))
			Tt=max(Tt);  	
	            if(is.real(Tt) && is.null(Tt)){  
		      # Upwind condition check, current distance must be larger
		      # then direct neighbours used in solution
		      DirectNeigbInSol=Tm[is.finite(Tm)]; #browser()
		      #if(nnz(DirectNeigbInSol>=Tt)>0) #if there's some element != zero
		            if( length(DirectNeigbInSol[DirectNeigbInSol>=Tt])>0 )
			                Tt=min(DirectNeigbInSol)+(1/(F[i,j]+eps));	            
		            }else{
			    		Tt=min(Tm)+(1/(F[i,j]+eps)); #when the root was complex 
		    }
	            # Update distance in trial list or add to neigbour list
	            if(Status[i,j]>0) # if it's in the trial list
	                Ltrial[1,Status[i,j]]=min(Tt,Ltrial[1,Status[i,j]])
	            else{
	                nLtrial=nLtrial+1;
			Status[i,j]=nLtrial;
	                # If running out of memory at a new block
	                if(nLtrial>neigh_free){
				neigh_free = neigh_free +100000; Ltrial[1,neigh_free]=0; }
	                Ltrial[,nLtrial]=c(Tt,i,j,timer);
	                T[i,j]=Tt;
	            }
	        }#if not Alive
	    }#for neighbors
	}#if timer
}
return (T)
}
