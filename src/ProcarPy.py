#-*- coding : utf-8 -*-
import numpy as np
import pylab as p
from matplotlib.lines import Line2D


class PROCARBandStructure:
	"""Class to parse and plot PROCAR VASP Data File"""
	def __init__(self, filename, path=None, ef=None, pathstyle="discontinuous",SO=False,spin=False):
		self.BandFile = filename
		self.label = []
		if path == None: self.Kpath = []
		else: self.Kpath = path
		self.Ef = ef
		self.SO = SO
		self.Spin = spin
		f = open(self.BandFile, "r")
		self.BandInfos = f.readlines()
		f.close()
		self.nkpts,self.nband,self.nions=[int(i) for i in self.BandInfos[1].split() if i.isdigit()]
		self.kpts = np.zeros((self.nkpts,3))
		if self.Spin: self.ene = np.zeros((self.nkpts,self.nband,2))
		else : self.ene = np.zeros((self.nkpts,self.nband,1))
		self.norbitals=len(self.BandInfos[7].rstrip("\n").split())-1

		index=3
		# Parse DATA for Spin-Orbit Coupling calculation
		if SO == True:
			self.ene_persentage = np.zeros((self.nkpts,self.nband,self.nions,4,self.norbitals))
			for ik in range(self.nkpts):
				if self.nions == 1: ind_k = index + ik*((4 + 4*(self.nions))*self.nband + 3)
				else: ind_k = index + ik*((4 + 4*(self.nions+1))*self.nband + 3)
				# self.kpts[ik] = np.array(self.BandInfos[ind_k].rstrip("\n").split()[3:6],float)
				k_info=self.BandInfos[ind_k][19:51] # to avoid negative sign "-0.00000000-0.00000000-0.00000000"
				self.kpts[ik] = np.array([k_info[i*11:i*11+11] for i in range(3)],float)
				for iband in range(self.nband):
					if self.nions == 1: ind_band = ind_k + 2 + iband*(2 + (4*(self.nions)+2))
					else: ind_band = ind_k + 2 + iband*(2 + (4*(self.nions+1)+2))
					self.ene[ik,iband,0] = float(self.BandInfos[ind_band].rstrip("\n").split()[4])
					for mx in range(4):
						if self.nions == 1:  ind_ion = ind_band + 3 + mx*(self.nions)
						else: ind_ion = ind_band + 3 + mx*(self.nions+1)
						for ion in range(self.nions):
							self.ene_persentage[ik,iband,ion,mx]=np.array(self.BandInfos[ind_ion+ion].rstrip("\n").split()[1:self.norbitals+1],float)
		#Parse DATA for no-Spin Calculation					
		elif self.Spin == False:
			self.ene_persentage = np.zeros((self.nkpts,self.nband,self.nions,self.norbitals))
			for ik in range(self.nkpts):
				if self.nions == 1: ind_k = index + ik*((4+(self.nions))*self.nband+3)
				else: ind_k = index + ik*((4+(self.nions+1))*self.nband+3)
				# self.kpts[ik] = np.array(self.BandInfos[ind_k].rstrip("\n").split()[3:6],float)
				k_info=self.BandInfos[ind_k][19:51] # to avoid negative sign "-0.00000000-0.00000000-0.00000000"
				self.kpts[ik] = np.array([k_info[i*11:i*11+11] for i in range(3)],float)
				for iband in range(self.nband):
					if self.nions == 1: 
						ind_band = ind_k + 2 + iband*(2 + (self.nions)+2)
					else:
						ind_band = ind_k + 2 + iband*(2 + (self.nions+1)+2)
					self.ene[ik,iband,0] = float(self.BandInfos[ind_band].rstrip("\n").split()[4])
					ind_ion = ind_band + 3
					for ion in range(self.nions):
						self.ene_persentage[ik,iband,ion]=abs(np.array(self.BandInfos[ind_ion+ion].rstrip("\n").split()[1:self.norbitals+1],float))
		#Parse DATA for Spin-Polarized calculations
		elif self.Spin == True:
			self.ene_persentage = np.zeros((self.nkpts,self.nband,self.nions,self.norbitals*2))
			if self.nions == 1: index_k_down = index + self.nkpts*((4+(self.nions))*self.nband+3) + 1
			else: index_k_down = index + self.nkpts*((4+(self.nions+1))*self.nband+3) + 1
			for ik in range(self.nkpts):
				if self.nions == 1:
					ind_k_up = index + ik*((4+(self.nions))*self.nband+3)
					ind_k_down = index_k_down + ik*((4+(self.nions))*self.nband+3)
				else:
					ind_k_up = index + ik*((4+(self.nions+1))*self.nband+3)
					ind_k_down = index_k_down + ik*((4+(self.nions+1))*self.nband+3)
				# self.kpts[ik] = np.array(self.BandInfos[ind_k_up].rstrip("\n").split()[3:6],float)
				k_info=self.BandInfos[ind_k_up][19:51] # to avoid negative sign "-0.00000000-0.00000000-0.00000000"
				self.kpts[ik] = np.array([k_info[i*11:i*11+11] for i in range(3)],float)
				for iband in range(self.nband):
					if self.nions == 1:
						ind_band_up = ind_k_up + 2 + iband*(2 + (self.nions)+2)
						ind_band_down = ind_k_down + 2 + iband*(2 + (self.nions)+2)
					else:
						ind_band_up = ind_k_up + 2 + iband*(2 + (self.nions+1)+2)
						ind_band_down = ind_k_down + 2 + iband*(2 + (self.nions+1)+2)
					self.ene[ik,iband,0] = float(self.BandInfos[ind_band_up].rstrip("\n").split()[4])
					self.ene[ik,iband,1] = float(self.BandInfos[ind_band_down].rstrip("\n").split()[4])
					ind_ion_up = ind_band_up + 3
					ind_ion_down = ind_band_down + 3
					for ion in range(self.nions):
						self.ene_persentage[ik,iband,ion,0:self.norbitals]=abs(np.array(self.BandInfos[ind_ion_up+ion].rstrip("\n").split()[1:self.norbitals+1],float))
						self.ene_persentage[ik,iband,ion,self.norbitals:2*self.norbitals+1]=abs(np.array(self.BandInfos[ind_ion_down+ion].rstrip("\n").split()[1:self.norbitals+1],float))

		del self.BandInfos

		if pathstyle == "continuous":
			self.WaveKpoints = np.zeros(self.nkpts)
			diff = np.diff(self.kpts, axis=0)
			norm = np.array([np.linalg.norm(diff[i]) for i in range(len(diff))])
			for i in range(1, self.nkpts):
				self.WaveKpoints[i] = self.WaveKpoints[i - 1] + norm[i - 1]
		else:
			self.WaveKpoints = np.arange(self.nkpts)

		if len(self.Kpath) > 0:
			self.Ksym = [int(self.nkpts / (len(path) - 1)) * i - 1 for i in range(1, len(path))]
			self.Ksym.insert(0, 0)
		else:
			self.Ksym = np.arange(self.nkpts)[::int(self.nkpts/6)]

	@staticmethod
	def OrbitalIndex(x,n,spin):
		if spin :
			if n == 4:
				Orbitals = {
					"s_up"   : 0,   "s_down"   : n ,
					"p_up"   : 1,   "p_down"   : n+1,
					"d_up"   : 2,   "d_down"   : n+2,
					"tot_up" : n-1, "tot_down" : 2*n-1,
					"up"     : 0,   "down"     : 1,
					}
			elif n == 10:
				Orbitals = {
					"s_up"      : 0,   "s_down"      : n ,
					"px_up"     : 3,   "px_down"     : n+3,
					"py_up"     : 1,   "py_down"     : n+1,
					"pz_up"     : 2,   "pz_down"     : n+2,
					"dxy_up"    : 4,   "dxy_down"    : n+4,
					"dyz_up"    : 5,   "dyz_down"    : n+5,
					"dz2_up"    : 6,   "dz2_down"    : n+6,
					"dxz_up"    : 7,   "dxz_down"    : n+7,
					"dx2-y2_up" : 8,   "dx2-y2_down" : n+8,
					"tot_up"    : n-1, "tot_down"    : 2*n-1,
					"up"        : 0,   "down"        : 1,
					}
		else:
			if n == 4:
				Orbitals = {
					"s"   : 0,
					"p"   : 1,
					"d"   : 2,
					"tot" : n-1,
					}
			elif n == 10:
				Orbitals = {
					"s"   : 0,
					"px"  : 3,   "py"  : 1, "pz"  : 2,
					"dxy" : 4,   "dyz" : 5, "dz2" : 6, "dxz" : 7, "dx2-y2"  : 8,
					"tot" : n-1,
					}
		return Orbitals[x]

	def getOrbitalBandsPlot(self, orbital=None, atom=1, magn="mtot", sign="+", marker="o", color="b", alpha=1, scale=100, label=None):
		'''
		This method allows to plot the Projected Band structure of atomic orbital (s, p, p_x, p_y, ...) for a defined atom.
		'''
		if orbital == None:
			if self.Spin : orbital = "tot_up"
			else : orbital = "tot"
		ind_orbital = self.OrbitalIndex(orbital,self.norbitals,self.Spin)
		SS = np.copy(self.ene_persentage)*scale
		atom -= 1

		if self.SO == True:
			magn_dict = {"mtot": 0, "mx": 1, "my": 2, "mz": 3}
			ind_magn = magn_dict[magn]
			if magn == "mtot":
				[self.ax.scatter(self.WaveKpoints, self.ene[:, i, 0] - self.Ef if self.Ef != None else self.ene[:, i, 0],alpha=alpha,marker=marker,c=color,s=SS[:,i,atom,ind_magn,ind_orbital]) for i in range(self.nband)]
				self.label.append(Line2D([0],[0],linestyle="",marker=marker,alpha=alpha,c=color,ms=10,label=label if label !=None else "Atom#"+str(atom+1)+" ("+magn+","+orbital+")"))
			else:
				percentage = SS[:,:,atom,ind_magn,ind_orbital]
				if sign == "+":
					posMagn = np.where(percentage > 0)
					self.ax.scatter(self.WaveKpoints[posMagn[0]],self.ene[posMagn[0],posMagn[1]] - self.Ef if self.Ef != None else self.ene[posMagn[0],posMagn[1]],marker=marker,color=color,s=percentage[posMagn])
					self.label.append(Line2D([0],[0],linestyle="",marker=marker,alpha=alpha,c=color,ms=10,label=label if label !=None else "Atom#"+str(atom+1)+" ("+magn+","+orbital+" $+$)"))
				elif sign == "-":
					negMagn = np.where(percentage < 0)	
					self.ax.scatter(self.WaveKpoints[negMagn[0]],self.ene[negMagn[0],negMagn[1]] - self.Ef if self.Ef != None else self.ene[negMagn[0],negMagn[1]],marker=marker,color=color,s=abs(percentage[negMagn]))
					self.label.append(Line2D([0],[0],linestyle="",marker=marker,alpha=alpha,c=color,ms=10,label=label if label !=None else "Atom#"+str(atom+1)+" ("+magn+","+orbital+" $-$)"))

		else:
			[self.ax.scatter(self.WaveKpoints, self.ene[:, i, 0] - self.Ef if self.Ef != None else self.ene[:, i, 0],marker=marker,alpha=alpha,c=color,s=SS[:,i,atom,ind_orbital]) for i in range(self.nband)]
			self.label.append(Line2D([0],[0],linestyle="",marker=marker,alpha=alpha,c=color,ms=10,label=label if label !=None else "Atom#"+str(atom+1)+" ("+orbital+")"))
	
	def getTotalBandsPlot(self, bandspin="up", lw=1, color="b", alpha=1, label=None):
		'''
		This method allows to plot the total electronic band structure from PROCAR file.
		'''
		if self.Spin :
			ind_spin = self.OrbitalIndex(bandspin,self.norbitals,self.Spin)
			if label == None: band_label = "Spin "+bandspin
			else: band_label = label
		else:
			ind_spin = 0
			if label == None: band_label = "Total Bands"
			else: band_label = label
		self.ax.plot(self.WaveKpoints, self.ene[:, :, ind_spin] - self.Ef if self.Ef != None else self.ene[:, :, ind_spin], c=color, alpha= alpha, lw=lw)
		self.label.append(Line2D([0],[0],c=color,alpha=alpha,ms=10,label=band_label))

	def getAtomsRangeBandsPlot(self, orbital=None, AtomRange=None, magn="mtot", sign="+", marker="o", color="b", alpha=1, scale=100, label=None):
		'''
		This method allows to plot the Projected Band structure of atomic orbital (s, p, p_x, p_y, ...) for a defined range of atoms.
		'''
		if AtomRange == None: AtomRange = np.arange(self.nions)
		else: AtomRange = np.array(AtomRange) - 1
		if orbital == None:
			if self.Spin : orbital = "tot_up"
			else : orbital = "tot"
		ind_orbital = self.OrbitalIndex(orbital,self.norbitals,self.Spin)
		if self.SO:
			magn_dict = {"mtot": 0, "mx": 1, "my": 2, "mz": 3}
			ind_magn = magn_dict[magn]
			Atoms_percentage = np.array([np.sum([self.ene_persentage[i,j,AtomRange,ind_magn,ind_orbital] for j in range(self.nband)],axis=1) * scale for i in range(self.nkpts)])
			if magn == "mtot":
				[self.ax.scatter(self.WaveKpoints, self.ene[:, i, 0] - self.Ef if self.Ef != None else self.ene[:, i, 0],marker=marker,alpha=alpha,c=color,s=Atoms_percentage[:,i]) for i in range(self.nband)]
				self.label.append(Line2D([0],[0],linestyle="",marker=marker,alpha=alpha,c=color,ms=10,label=label if label !=None else "Atom#"+str(np.array(AtomRange)+1)+" ("+magn+","+orbital+")"))
			else:
				if sign == "+":
					posMagn = np.where(Atoms_percentage > 0)
					self.ax.scatter(self.WaveKpoints[posMagn[0]],self.ene[posMagn[0],posMagn[1]] - self.Ef if self.Ef != None else self.ene[posMagn[0],posMagn[1]],marker=marker,color=color,s=Atoms_percentage[posMagn])
					self.label.append(Line2D([0],[0],linestyle="",marker=marker,alpha=alpha,c=color,ms=10,label=label if label !=None else "Atom#"+str(np.array(AtomRange)+1)+" ("+magn+","+orbital+" $+$)"))
				elif sign == "-":
					negMagn = np.where(Atoms_percentage < 0)	
					self.ax.scatter(self.WaveKpoints[negMagn[0]],self.ene[negMagn[0],negMagn[1]] - self.Ef if self.Ef != None else self.ene[negMagn[0],negMagn[1]],marker=marker,color=color,s=abs(Atoms_percentage[negMagn]))
					self.label.append(Line2D([0],[0],linestyle="",marker=marker,alpha=alpha,c=color,ms=10,label=label if label !=None else "Atom#"+str(np.array(AtomRange)+1)+" ("+magn+","+orbital+" $-$)"))

		else:
			Atoms_percentage = np.array([np.sum([self.ene_persentage[i,j,AtomRange,ind_orbital] for j in range(self.nband)],axis=1) * scale for i in range(self.nkpts)])
			[self.ax.scatter(self.WaveKpoints, self.ene[:, i, 0] - self.Ef if self.Ef != None else self.ene[:, i, 0],alpha=alpha,marker=marker,c=color,s=Atoms_percentage[:,i]) for i in range(self.nband)]
			self.label.append(Line2D([0],[0],linestyle="",marker=marker,alpha=alpha,c=color,ms=10,label=label if label !=None else "Atom#"+str(np.array(AtomRange)+1)+" ("+orbital+")"))	
	
	@property
	def getKpoints(self):
		return self.kpts

	@property
	def getEnergies(self):
		return self.ene

	def getBandsGap(self, band1, band2, bandspin=""):
		'''
		This method allows to calculate the band gap between two defined bands.
		def getBandsGap(self,band1,band2,bandspin):
		- band1: sets the first band index.
		- band2: sets the second band index.
		- bandspin: sets the band spin ("up", "down") if spin = True.
		'''
		bandspin_index = {"": 0, "up": 0, "down": 1}
		ind_s = bandspin_index[bandspin]
		E1 = np.amax(self.ene[:,band1-1,ind_s])
		E2 = np.amin(self.ene[:,band2-1,ind_s])
		gap = max(E2-E1,0.0)
		print("Band Gap between band# {} & band# {} = {} eV".format(band1,band2,gap))
		return gap

	def Init_Fig(self,width=8,height=6):
		self.label = []
		p.close()
		f, self.ax = p.subplots(figsize=(width,height))
		f.subplots_adjust(right=0.75)

	def PlotShow(self, ymin=None, ymax=None, xmin=None, xmax=None, savefile=None):
		if ymin == None: ymin = self.ax.get_ylim()[0]
		if ymax == None: ymax = self.ax.get_ylim()[1]
		if xmin == None: xmin = self.WaveKpoints[0]
		if xmax == None: xmax = self.WaveKpoints[-1]
		self.ax.set_ylim(ymin, ymax)		
		self.ax.set_ylabel("$\epsilon - \epsilon_{Fermi}$ (eV)" if self.Ef != None else "$Energy$ (eV)")
		self.ax.set_xlabel("$k-$points")
		self.ax.set_xticks(self.WaveKpoints[self.Ksym])
		self.ax.set_xticklabels(self.Kpath if len(self.Kpath) > 0 else np.round(self.WaveKpoints[self.Ksym], 2))
		self.ax.set_xlim(xmin, xmax)

		# self.ax.legend(handles=self.label,loc="best")
		self.ax.legend(handles=self.label,bbox_to_anchor=(1.04,1), loc="upper left")
		if self.Ef != None:
			self.ax.hlines(0, self.WaveKpoints[0], self.WaveKpoints[-1], linestyle="dashed")
		if len(self.Kpath) > 0:
			self.ax.vlines(self.WaveKpoints[self.Ksym], ymin, ymax, "gray", alpha=0.75)
		if savefile == None:
			p.show()
		else:
			p.savefig(savefile)	

	def __repr__(self):
		print("Bands File : {}".format(self.BandFile))
		print("nkpts : {}, nband : {} , nions : {} ".format(self.nkpts, self.nband, self.nions))
		print("Spin calcuations : {}".format(self.Spin))
		print("Spin Orbit calcuations : {}".format(self.SO))
		print("K-Path : {}".format(self.Kpath))
		print("Efermi : {}".format(self.Ef))
		return "PROCARBandStructure('{}')".format(self.BandFile)

	def __str__(self):
		return self.__repr__()