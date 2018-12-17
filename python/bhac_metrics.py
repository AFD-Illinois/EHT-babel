
import numpy as np

def gcov_logks_2d(x1,x2,a):
  """ Returns logks g_munu(x1,x2) for given a.
      Courtesy of (stolen from) George Wong"""
  r = np.exp(x1)
  th = x2
  r,th = [ x.T for x in np.meshgrid(r,th) ]
  
  # calculate normal ks
  gcov = np.zeros((len(x1),len(x2),4,4))
  cth = np.cos(th)
  sth = np.sin(th)
  s2 = sth*sth
  rho2 = r*r + a*a*cth*cth
  gcov[:,:,0,0] = (-1. + 2. * r / rho2)
  gcov[:,:,0,1] = (2. * r / rho2)
  gcov[:,:,0,3] = (-2. * a * r * s2 / rho2)
  gcov[:,:,1,0] = gcov[:,:,0,1]
  gcov[:,:,1,1] =  (1. + 2. * r / rho2)
  gcov[:,:,1,3] =  (-a * s2 * (1. + 2. * r / rho2))
  gcov[:,:,2,2] = rho2
  gcov[:,:,3,0] = gcov[:,:,0,3]
  gcov[:,:,3,1] = gcov[:,:,1,3]
  gcov[:,:,3,3] = s2 * (rho2 + a*a * s2 * (1. + 2. * r / rho2))
  # transform to logks
  trans = np.zeros([len(x1),len(x2),4,4])
  trans[:,:,0,0] = 1.
  trans[:,:,1,1] = r
  trans[:,:,2,2] = 1.
  trans[:,:,3,3] = 1.
  #return np.transpose(np.matmul(trans,np.matmul(gcov,trans)),(1,0,2,3))
  return np.matmul(trans,np.matmul(gcov,trans))


