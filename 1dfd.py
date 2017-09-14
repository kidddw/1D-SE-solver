#!/usr/bin/python -tt

#add other potentials


"""One-dimensional Schrodinger equation solver

   choose basis 
"""

import sys
import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt


class Define_Parameters:
  def __init__(self, unit):
    if unit != "atomic" and unit != "eV": 
      print "must be either (1) atomic or (2) eV"
      sys.exit()
    if unit == "eV":
      self.units = unit
      self.hbar = 0.658211899 # Planck constant over 2*pi in eV*fs
      self.CoulombConst = 14.39964415 # 1/(4*pi*eps_0), where eps_0=0.005526349868 eV/(Angstrom*V^2)
      self.h2m = 3.80998174348 # (hbar^2)/(2*m_e) in eV*Angstom^2
    if unit == "atomic":
      self.units = unit
      self.hbar = 1.0
      self.CoulombConst = 1.0
      self.h2m = 0.5
    self.emass = 0.5*self.hbar*self.hbar/self.h2m # electron mass


class Real_Space_Grid:
  def __init__(self, a, b, n, param):
    """Define Hamiltonian and interpret eigensolutions for a real-space grid representation
    9-point finite difference is employed"""
    print "REAL-SPACE GRID BASIS"
    self.a = a
    self.b = b
    self.n = n
    self.grid_spacing = float(b - a)/(n - 1.0)
    self.length = (a - b)
    self.grid_point = []
    for i in range(n):
      x = i*self.grid_spacing + a
      self.grid_point.append(x)
    self.H_mat = self.define_fd9(param) + self.potential(param,"harmonic oscillator")
    self.O_mat = np.identity(n) # overlap is identity


  def define_fd9(self,param):
    """Initialize a kinetic matrix using 9-point finite difference"""
#    https://en.wikipedia.org/wiki/Finite_difference_coefficient
    fd_coeff = []
    fd_coeff.append(-205.0/72.0) 
    fd_coeff.append(8.0/5.0)
    fd_coeff.append(-1.0/5.0)
    fd_coeff.append(8.0/315.0)
    fd_coeff.append(-1.0/560.0)
    arr =[]
    for i in range(self.n):
      row= []
      for j in range(self.n):
        for k in range(5):
          if j==i+k or j==i-k: row.append(fd_coeff[k])
        if abs(j-i) > 4: row.append(0.0)
      arr.append(row)
    mat = np.matrix(arr)
    mat = -1.0 * param.h2m / self.grid_spacing**2 * mat
    return mat


  # choices of harmonic oscillator, ...
  # To Do: 
  #   - add more choices of potential
  #   - add read in of omega for harmonic osc.
  def potential(self,param,form):
    """Define a potential matrix"""
    if form != "harmonic oscillator" and form != "morse": 
      print "must be either (1) harmonic oscillator or (2) morse"
      sys.exit()
    pot = []
    if form == "harmonic oscillator": pot = define_ho(self,param,1.0)
#    if form == "(...)": pot = define_(...)
    I = np.identity(self.n)
    mat = np.diag(pot)
    return mat


  def normalize(self,psi):
    for ind in range(psi.shape[1]):
      su = 0.0
      for i in range(psi.shape[0]): 
        su = su + psi[i,ind] * psi[i,ind] * self.grid_spacing
      print su
      for i in range(psi.shape[0]):
        psi[i,ind] = psi[i,ind] / np.sqrt(su)
    return psi


  # write out first energies and plot ground state wavefunction
  def output_eigensolutions(self,en,wf,n,param):
    """Interpret eigensolutions"""
    for item in en[0:n]:
      print item
    i = 0
    while (i < n):
      plt.plot(self.grid_point,wf[:,i])
      i = i + 1
    plt.ylabel('$\psi (x)$')
    if param.units == "eV": plt.xlabel('$x$ [Angstrom]')
    if param.units == "atomic": plt.xlabel('$x$ [a.u.]')
    plt.draw()



class Gaussian:
  def __init__(self, x0, a0, n, param):
    """Define Hamiltonian and interpret eigensolutions for a gaussian basis representation
    Gaussians are all centered at x=0"""
    print "VARIABLE WIDTH GAUSSIAN BASIS"
    self.x0 = x0
    self.a0 = a0
    self.n = n
    self.nu = []
    for i in range(n):
      x = 1.0 / (a0 * x0**i)**2
      self.nu.append(x)
    self.O_mat = self.Overlap_mat()
    self.H_mat = self.T_mat(param) + self.HO_mat(param,1.0)


  def Overlap_mat(self):
    """Overlap matrix for the Gaussian basis where all Gaussians are centered at x=0"""
    arr =[]
    for i in range(self.n):
      row= []
      for j in range(self.n):
        num = 2.0 * np.sqrt(self.nu[i] * self.nu[j])
        den = self.nu[i] + self.nu[j]
        mat_elem = np.sqrt( num / den )
        row.append(mat_elem)
      arr.append(row)
    mat = np.matrix(arr)
    return mat


  def T_mat(self,param):
    """Kinetic matrix for the Gaussian basis where all Gaussians are centered at x=0"""
    arr = []
    for i in range(self.n):
      row= []
      for j in range(self.n):
        mat_elem = 2.0 * (self.nu[i] * self.nu[j]) / (self.nu[i] + self.nu[j])
        mat_elem = mat_elem * self.O_mat[i,j]
        row.append(mat_elem)
      arr.append(row)
    mat = param.h2m * np.matrix(arr)
    return mat


  def HO_mat(self,param,omega):
    """Harmonic oscillator potential matrix for the Gaussian basis where 
    all Gaussians are centered at x=0"""
    fac = 0.5 * param.emass * omega**2
    arr = []
    for i in range(self.n):
      row= []
      for j in range(self.n):
        mat_elem = fac / (2.0 * (self.nu[i] + self.nu[j])) * self.O_mat[i,j]
        row.append(mat_elem)
      arr.append(row)
    mat = np.matrix(arr)
    return mat


  def normalize(self,psi):
    for ind in range(self.n):
      su = 0.0
      for i in range(self.n):
        for j in range(self.n): 
          su = su + psi[i,ind] * psi[j,ind] * self.O_mat[i,j]
      for i in range(self.n):
        psi[i,ind] = psi[i,ind] / np.sqrt(su)
    return psi


  # write out first energies and plot ground state wavefunction
  def output_eigensolutions(self, en, wf, n, param):
    """Interpret eigensolutions. Hardcoded in -5 to 5"""
    for item in en[0:n]:
      print item
    a = -5.0
    b = 5.0
    n_points = 100
    h = float(b - a)/(n_points - 1.0)
    grid_point = []
    for ind in range(n_points):
      x = ind*h + a
      grid_point.append(x)
    i = 0
    while (i < n):
      phi = []
      for ind in range(n_points):
        su = 0.0
        for j in range(self.n):
          su = su + wf[j,i] * self.gauss_func(self.nu[j], grid_point[ind])
        phi.append(su)
      plt.plot(grid_point,phi)
      i = i + 1
    plt.ylabel('$\psi (x)$')
    if param.units == "eV": plt.xlabel('$x$ [Angstrom]')
    if param.units == "atomic": plt.xlabel('$x$ [a.u.]')
    plt.draw()

       
  def gauss_func(self,nu,x):
    """Gaussian basis function in real-space"""
    norm_const = np.sqrt(np.sqrt(2.0 * nu / np.pi))
    gauss = norm_const * np.exp(-1.0 * nu * x**2)
    return gauss


# Harmonic Oscillator potential
# https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator
def define_ho(grid,param,omega): 
  ho_arr = []
  for i in range(grid.n): 
    val = 0.5 * param.emass * omega**2 * grid.grid_point[i]**2 
    ho_arr.append(val)
  return ho_arr



# Initializes hamiltonian, solves eigenproblem, outputs first three solutions
def main():

  # Declare units
  myparam = Define_Parameters("atomic")

  # Use real-space grid basis
#  basis = Real_Space_Grid(-5, 5, 100, myparam)

  # Use variable-width Gaussian basis
  basis = Gaussian(1.14, 0.01, 100, myparam)

  # Solve generalized eigenvalue problem
  energy, psi = eigh(basis.H_mat, basis.O_mat)

  # Output solutions
  psi = basis.normalize(psi)
  basis.output_eigensolutions(energy, psi, 3, myparam) 

  plt.show()

if __name__ == '__main__':
  main()




