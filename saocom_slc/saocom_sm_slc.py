from osgeo import gdal
import os
import rasterio
from xml.etree import ElementTree as ET
import numpy as np
from datetime import datetime,timedelta
import math
import subprocess
import shutil
import glob
import warnings
warnings.filterwarnings('ignore')

###################
'''
This script uses a class saocom_sm_slc and a function read_saocom that are used to read SAOCOM-1 Stripmap Data in SLC format
to GMTSAR.  To test it, the script must be located with in the same folder of the input data (xemt files and unzipped data folders).
'''
###################

def read_saocom(polarizations = None):
	for filename in os.listdir('.'):
		if filename.endswith('.xemt'):
			a = saocom_sm_slc(filename)
			a.process_xml(polarizations)
			
class saocom_sm_slc:
	def __init__(self, xemt_file):
		self.xemt = xemt_file
		
	def get_root(self, xfile):
		f = open(xfile,'rt')
		tree =  ET.parse(f)
		root =  tree.getroot()
		f.close()
		return root
		
	def  get_data_files(self,polarizations):
		component_dir = self.xemt.split('.')[0]
		root = self.get_root(self.xemt)
		product = root.findall('product')
		datafile =  [feat.findall("dataFile") for feat in product][0]
		components = [feat.findall("components") for feat in datafile][0]
		component = [os.path.join(component_dir,feat[0].text) for feat in components[0] if feat[0].text.startswith("Data/slc-acq")]
		
		#Filter the list by the polarizations selected by the user
		if polarizations:
			pol_list = [pol.lower() for pol in polarizations] # lowercase all the letters in case the user's input is in uppercase
			component = [c for c in component if c.split('.')[0][-2:] in pol_list]
		return component
		
	def process_xml(self, polarizations):
		xml_list = self.get_data_files(polarizations)
		for xml_file in xml_list:
			slc_file = xml_file.split('.')[0]
			print(f'processing file {slc_file}') 
			self.read_params_xml(xml_file)
			self.write_led()
			self.write_prm()
			self.write_slc(slc_file)
			os.system('chmod 777 -R *')
			self.move_files(self.pol)
			
	def read_params_xml(self, xml_file):
		self.root = self.get_root(xml_file)
		self.Channel = self.root.findall('Channel')
		self.RasterInfo = [feat.findall("RasterInfo") for feat in self.Channel][0]
		self.DataSetInfo = [feat.findall("DataSetInfo") for feat in self.Channel][0]
		self.SwathInfo =  [feat.findall("SwathInfo") for feat in self.Channel][0]
		self.SamplingConstants =  [feat.findall("SamplingConstants") for feat in self.Channel][0]
		self.BurstInfo =  [feat.findall("BurstInfo") for feat in self.Channel][0]
		self.StateVectorData = [feat.findall("StateVectorData") for feat in self.Channel][0]
		self.DopplerCentroid = [feat.findall("DopplerCentroid") for feat in self.Channel][0]
		self.Pulse = [feat.findall("Pulse") for feat in self.Channel][0]
		self.Burst = [feat.findall("Burst") for feat in self.BurstInfo][0] # For SAOCOM Stripmap, there is only 1 Burst

		#Get the constant parameters for the PRM file
		self.prm_constants()
		
		#Get the input file basename
		self.input_file = xml_file.split('.')[0]
		
		#Read the parameters from the XML file:
		self.SensorName = [feat.find("SensorName").text for feat in self.DataSetInfo][0]
		self.Polarization = [feat.find("Polarization").text for feat in self.SwathInfo][0]
		self.rng_sampling_rate = [feat.find("frg_hz").text for feat in self.SamplingConstants][0]
		self.fc_hz = float([feat.find("fc_hz") for feat in self.DataSetInfo][0].text)
		self.pulse_dur = float([feat.find("PulseLength").text for feat in self.Pulse][0])
		self.chirp_slope = float([feat.find("Bandwidth").text for feat in self.Pulse][0]) / self.pulse_dur
		self.prf = float([feat.find("faz_hz").text for feat in self.SamplingConstants][0])
		self.l = self.c / self.fc_hz
		self.SamplesStart = float([feat.find("SamplesStart").text for feat in self.RasterInfo][0])
		self.near_range = self.c * self.SamplesStart / 2.
		self.AzimuthStartTime = [feat.find("AzimuthStartTime").text for feat in self.Burst][0]
		self.dateTime,self.decSeconds = self.AzimuthStartTime.split('.')
		self.microsec = float("0."+self.decSeconds)*1e6
		self.dt = datetime.strptime(self.dateTime,'%d-%b-%Y %H:%M:%S')
		self.dt = self.dt + timedelta(microseconds=self.microsec)
		self.hour = self.dt.hour
		self.minute = self.dt.minute
		self.second = self.dt.second
		self.micro = self.dt.microsecond
		self.clock_start = ((self.hour/24)+(self.minute/1440)+(self.second/86400) + self.micro/(1000000*86400)) + int(self.dt.strftime('%j'))
		self.SC_clock_start = self.clock_start + self.dt.year * 1000
		self.day = int(self.dt.strftime('%j'))
		self.orbdir = [feat.find("OrbitDirection").text for feat in self.StateVectorData][0]
		self.lookdir = [feat.find("SideLooking").text for feat in self.DataSetInfo][0] 
		self.num_lines =  int([feat.find("Lines").text for feat in self.RasterInfo][0])
		self.nrows = self.num_lines
		self.num_valid_az = self.num_lines
		self.num_rng_bins = int([feat.find("Samples").text for feat in self.RasterInfo][0])
		self.bytes_per_line = self.num_rng_bins * 4
		self.good_bytes = self.bytes_per_line
		self.LinesStep = float([feat.find("LinesStep").text for feat in self.RasterInfo][0]) # LinesStep: Used to compute clock_stop
		self.clock_stop = self.clock_start + self.num_lines * self.LinesStep / 86400
		self.SC_clock_stop = self.SC_clock_start + self.num_lines * self.LinesStep / 86400
		
		#Use Polarization , SensorName and DateTime to form the file basename
		self.pol = self.Polarization.split('/')[0]+self.Polarization.split('/')[1]
		self.base = (self.SensorName + "_" + str(self.dt.year) + str(self.dt.month).zfill(2) + str(self.dt.day).zfill(2)) + "_" + self.pol
		
		#Doppler Centroid (fd1 and fdd1 in PRM)
		
		dr = 0.5*self.c/float(self.rng_sampling_rate) #Theoretic range spacing (m)
		sample_step = float([feat.find("SamplesStep").text for feat in self.RasterInfo][0])
		
	# Get t0 in range direction
		trg_0 = []
		for feat in self.DopplerCentroid:
			for feat2 in feat.findall("trg0_s"):
				trg_0.append(float(feat2.text))
		
		# Get the t0 en azimuth direction
		taz_0 = []
		for feat in self.DopplerCentroid:
			for feat2 in feat.findall("taz0_Utc"):
				date = datetime.strptime(feat2.text.split('.')[0],'%d-%b-%Y %H:%M:%S')
				date_num = date.hour/24+date.minute/1440+date.second/86400 + int(date.strftime('%j')) + date.year * 1000
				taz_0.append(date_num)
				
		#Get the DC polynomials	
		dopRngTime = []
		for feat in self.DopplerCentroid:
			for feat2 in feat.findall("pol"):
				for val in feat2.findall("val"):
					dopRngTime.append(float(val.text))
					
		#Create an array in range direction based on the polynomials to get FDD1
		dc = []
		dc_median = []
		ddc = []
		ddc_median = []
		for i in range(len(trg_0)):
			dopRngTime_i = dopRngTime[7*i:7*i+7]
			dc_i = np.zeros(self.num_rng_bins)
			ddc_i = np.zeros(self.num_rng_bins)
			for j in range(self.num_rng_bins):
				trg = j*sample_step+trg_0[i]
				dc_i[j] = dopRngTime_i[0] + (trg-trg_0[i])*dopRngTime_i[1] + (taz_0[i] - taz_0[i])*dopRngTime_i[2] + (taz_0[i] - taz_0[i])*(trg-trg_0[i])*dopRngTime_i[3] + (trg-trg_0[i])**2*dopRngTime_i[4] + (trg-trg_0[i])**3*dopRngTime_i[5] + (trg-trg_0[i])**4*dopRngTime_i[6]
			dc.append(dc_i)
			dc_median_i = np.median(dc_i)
			dc_median.append(dc_median_i)
			for j in range(1,self.num_rng_bins):
				ddc_i[j-1] = dc_i[j]-dc_i[j-1]
			ddc.append(ddc_i)
			ddc_median_i = np.median(ddc_i)
			ddc_median.append(ddc_median_i)
		
		self.fd1 = np.mean(dc_median)
		self.fdd1 = np.mean(ddc_median) / dr
		
		#Read information for the LED file:
		
		self.nSV_n = int([feat.find("nSV_n") for feat in self.StateVectorData][0].text)
		self.year = self.dt.year
		
		#Get the t0 and duration of the acquisition of the image in seconds
		self.t0_image=math.trunc(self.hour*3600+self.minute*60+self.second+self.micro/1000000)
		self.dt_image = round(self.num_lines/self.prf)
		
		# Extract the stateVectors starting time and dtSV_s (separation in sec between each StateVector)
		self.t_ref_Utc =  [feat.find("t_ref_Utc") for feat in self.StateVectorData][0].text
		fmt = '%d-%b-%Y %H:%M:%S.%f000000'
		self.dt_orb = datetime.strptime(self.t_ref_Utc, fmt)
		self.year_orb = self.dt_orb.year
		self.day_orb = int(self.dt_orb.strftime('%j'))
		self.hour_orb = self.dt_orb.hour
		self.minute_orb = self.dt_orb.minute
		self.second_orb = self.dt_orb.second
		self.micro_orb = self.dt_orb.microsecond
		self.t0_orb=math.trunc(self.hour_orb*3600+self.minute_orb*60+self.second_orb+self.micro_orb/1000000)
		self.dtSV_s =  int([feat.find("dtSV_s") for feat in self.StateVectorData][0].text)
		
		#Create the acquisition and StateVectors time windows as arrays:
		self.t = [0]*self.dt_image
		self.t_orb = [0]*self.nSV_n
		
		for i in range(self.dt_image):
			self.t[i] = int(self.t0_image + self.dtSV_s*i)
		
		for i in range(self.nSV_n):
			self.t_orb[i] = int(self.t0_orb + self.dtSV_s * i)
			
		#Read the position and velocity StateVectors lists:
		self.pSV_m = 	[feat.find("pSV_m") for feat in self.StateVectorData][0]
		self.vSV_mOs = 	[feat.find("vSV_mOs") for feat in self.StateVectorData][0]
		
		
	def prm_constants(self):
		self.c = 299792458.0 
		self.nlooks = 1  # El numero de looks en SLC es 1
		self.rshift = 0
		self.ashift = 0
		self.sub_int_r = 0.0
		self.sub_int_a = 0.0
		self.stretch_r = 0.0
		self.stretch_a = 0.0
		self.a_stretch_r = 0.0
		self.a_stretch_a = 0.0
		self.first_sample = 1
		self.st_rng_bin = 1
		self.dtype = "a"
		self.SC_identity = 20 
		self.ra = 6378137.00  # equatorial radius
		self.rc = 6356752.31  # plar radius
		self.fddd1 = 0.0
		self.num_patches = 1
		self.chirp_ext = 0
		self.SLC_scale = 1
	
	def write_led(self):
		f = open(self.base + ".LED", "a")
		f.truncate(0)
		f.write('%s %s %s %.3f %.3f \n' % (self.nSV_n, self.year_orb, self.day_orb, self.t0_orb, self.dtSV_s))
		for i in range(self.nSV_n):
			s = self.t0_orb + self.dtSV_s * i
			x = round(float(self.pSV_m[i * 3].text), 6)
			y = round(float(self.pSV_m[i * 3 + 1].text), 6)
			z = round(float(self.pSV_m[i * 3 + 2].text), 6)
			vx = round(float(self.vSV_mOs[i * 3].text), 8)
			vy = round(float(self.vSV_mOs[i * 3 + 1].text), 8)
			vz = round(float(self.vSV_mOs[i * 3 + 2].text), 8)
			f.write('%s %s %.6f %.6f %.6f %.6f %.8f %.8f %.8f \n' % (self.year_orb, self.day_orb, s, x, y, z, vx, vy, vz))
		f.close()
	
	def write_prm(self):
		f = open(self.base + ".PRM", "a")
		f.truncate(0) 
		print("num_valid_az =  " + str(self.num_valid_az), file=f)
		print("nrows =  " + str(self.nrows), file=f)
		print("first_line =  " + str(1), file=f)
		print("deskew =  n", file=f)
		print("caltone =  0.0000", file=f)
		print("st_rng_bin =  1", file=f)
		print("flip_iq =  n", file=f)
		print("offset_video =  1", file=f)
		print("az_res =  0.0", file=f)
		print("nlooks =  1", file=f)
		print("chirp_ext =  0", file=f)
		print("scnd_rng_mig =  0", file=f)
		print("rng_spec_wgt =  1", file=f)
		print("rm_rng_band =  0.20", file=f)
		print("rshift =  " + str(self.rshift), file=f)
		print("ashift =  " + str(self.ashift), file=f)
		print("stretch_r =  " + str(self.stretch_r), file=f)
		print("stretch_a =  " + str(self.stretch_a), file=f)
		print("a_stretch_r =  " + str(self.a_stretch_r), file=f)
		print("a_stretch_a =  " + str(self.a_stretch_a), file=f)
		print("first_sample =  " + str(self.first_sample), file=f)
		print("SC_identity =  " + str(self.SC_identity), file=f)
		print("rng_samp_rate =  " + str(self.rng_sampling_rate), file=f)
		print("input_file =  " + self.input_file, file=f)
		print("num_rng_bins =  " + str(self.num_rng_bins), file=f)
		print("bytes_per_line =  " + str(self.bytes_per_line), file=f)
		print("good_bytes_per_line =  " + str(self.good_bytes), file=f)
		print("PRF =  " + str(self.prf), file=f)
		print("pulse_dur =  " + str(self.pulse_dur), file=f)
		print("near_range =  " + str(self.near_range), file=f)
		print("num_lines =  " + str(self.num_lines), file=f)
		print("num_patches =  " + str(self.num_patches), file=f)
		print("SC_clock_start =  " + str(self.SC_clock_start), file=f)
		print("SC_clock_stop =  " + str(self.SC_clock_stop), file=f)
		print("clock_start =  " + str(self.clock_start), file=f)
		print("clock_stop =  " + str(self.clock_stop), file=f)
		print("led_file =  " + self.base + ".LED", file=f)
		print("orbdir =  " + str(self.orbdir), file=f)
		print("lookdir =  " + str(self.lookdir), file=f)
		print("radar_wavelength =  " + str(self.l), file=f)
		print("chirp_slope =  " + str(self.chirp_slope), file=f)
		print("chirp_slope =  " + str(self.chirp_slope), file=f)
		print("I_mean =  0", file=f)
		print("Q_mean =  0", file=f)
		print("equatorial_radius =  " + str(self.ra), file=f)
		print("polar_radius =  " + str(self.rc), file=f)
		print("fd1 =  " + str(self.fd1), file=f)
		print("fdd1 =  " + str(self.fdd1), file=f)
		print("fddd1 =  " + str(self.fddd1), file=f)
		print("sub_int_r =  " + str(self.sub_int_r), file=f)
		print("sub_int_a =  " + str(self.sub_int_a), file=f)
		print("SLC_file =  " + self.base + ".SLC", file=f)
		print("dtype =  " + str(self.dtype), file=f)
		print("SLC_scale =  " + str(self.SLC_scale), file=f)
		
		f.close()
		
		#Write outputs from calc_dop_orb to the PRM (this will be required for coregistration)
		args = "calc_dop_orb " + self.base + ".PRM doppler.txt 0 " + str(self.fd1)
		os.system(args)
		os.system('grep SC_vel doppler.txt >> ' + self.base + ".PRM")
		os.system('grep earth_radius doppler.txt >> ' + self.base + ".PRM")
		os.system('grep SC_height doppler.txt >> ' + self.base + ".PRM")
		os.system("rm -f " + 'doppler.txt')
		os.system('chmod 777 -R *')
		

	def write_slc(self,slc_file):
		datafile = gdal.Open(slc_file)
		driver = datafile.GetDriver()
		cols = datafile.RasterXSize 
		rows = datafile.RasterYSize 
		data= datafile.GetRasterBand(1).ReadAsArray()
		datafile = None
		i = data.real*10000
		q = data.imag*10000
		outImag = driver.Create(self.base+"_iq.tif", cols, rows, 2, gdal.GDT_Int16)
		outBandImag = outImag.GetRasterBand(1)
		outDataImag = i
		outBandImag.WriteArray(outDataImag)
		outBandImag = outImag.GetRasterBand(2)
		outDataImag = q
		outBandImag.WriteArray(outDataImag)
		del outImag
		datafile = rasterio.open(self.base+'_iq.tif')
		data = datafile.read((1, 2))
		arr = np.transpose(data,[1,2,0])
		arr.tofile(self.base+'.SLC')
		del arr,i,q,data
		os.remove(self.base+'_iq.tif')
		
	def move_files(self,pol):
		if not os.path.exists(pol):
			os.mkdir(pol)
		else:
			for filename in os.listdir(pol):
				if filename.split('.')[0] == self.base:
					os.remove(os.path.join(pol,filename))
			
		SLC_files = glob.glob("*" + pol + ".SLC", recursive = True)
		PRM_files = glob.glob("*" + pol + ".PRM", recursive = True)
		LED_files = glob.glob("*" + pol + ".LED", recursive = True)
			
		for prm in PRM_files:
			shutil.move(os.path.join(os.getcwd(),prm),os.path.join(os.getcwd(),pol))
		for led in LED_files:
			shutil.move(os.path.join(os.getcwd(),led),os.path.join(os.getcwd(),pol))
		for slc in SLC_files:
			shutil.move(os.path.join(os.getcwd(),slc),os.path.join(os.getcwd(),pol))		
