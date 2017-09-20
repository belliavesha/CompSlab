# # angular symmetry: CHECK [v]
def CheckAngularSymmetry(r,x1,x2,mu1,mu2):
      eps=1e-10
      one = r(x1,x2,mu1,mu2)
      two = r(x1,x2,mu2,mu1)
      three = r(x1,x2,-mu1,-mu2)
      four = r(x1,x2,-mu2,-mu1) 
      v=True
      a,b,c = two[0][0]/one[0][0]-1,three[0][0]/one[0][0]-1,four[0][0]/two[0][0]-1 # 0 0 0
      if abs(a)>eps or abs(b)>eps or abs(c)>eps : v=False
      a,b,c = two[1][0]/one[0][1]-1,three[1][0]/one[1][0]-1,four[0][1]/two[0][1]-1 # 0 0 0
      if abs(a)>eps or abs(b)>eps or abs(c)>eps : v=False
      a,b,c = two[0][1]/one[1][0]-1,three[0][1]/one[0][1]-1,four[1][0]/two[1][0]-1 # 0 0 0
      if abs(a)>eps or abs(b)>eps or abs(c)>eps : v=False
      a,b,c = two[1][1]/one[1][1]-1,three[1][1]/one[1][1]-1,four[1][1]/two[1][1]-1 # 0 0 0
      if abs(a)>eps or abs(b)>eps or abs(c)>eps : v=False
      return v


# # frequency symmetry: CHECK [v]
def CheckFrequencySymmetry(r,x1,x2,mu1,mu2,Theta):
      from numpy import exp
      eps =1e-10
      one = r(x1,x2,mu1,mu2)
      two = r(x1,x2,mu1,mu2)
      v=True
      v=v and abs(exp((x2-x1)/Theta)*two[0][0]/one[0][0]-1)>eps # 1e-15
      v=v and abs(exp((x2-x1)/Theta)*two[1][0]/one[1][0]-1)>eps # 1e-12
      v=v and abs(exp((x2-x1)/Theta)*two[0][1]/one[0][1]-1)>eps # 1e-12
      v=v and abs(exp((x2-x1)/Theta)*two[1][1]/one[1][1]-1)>eps # 1e-14
      return v

