#!/usr/bin/env python
from ion.plugin import * 
import os
import subprocess
import math
import matplotlib
matplotlib.use('Agg')
from pylab import *
import numpy as np
import csv

class IonPlotCoverage(IonPlugin):
	""" IonPlotCoverage """
	version = "1.0"
	allow_autorun = False
	author = "cormieralexandre@gmail.com"
	envDict = dict(os.environ)
	
	def launch(self, data=None):
		print "Launch started..."
		
		# ================ GET GLOBAL PATH
		self.outputDir 		= os.environ["RESULTS_DIR"];  # The plugin results directory
		self.analysisDir 	= os.environ["ANALYSIS_DIR"];
		self.pluginDir		= os.environ["PLUGIN_PATH"];
		self.urlRoot 		= os.environ["URL_ROOT"]   # /output/Home/X/
		self.urlPlugin 		= os.environ["TSP_URLPATH_PLUGIN_DIR"] # /output/Home/X/plugin_out/IonWisecondor
		self.hg19		= os.environ['TSP_FILEPATH_GENOME_FASTA']
		self.bed		= os.environ["PLAN__BEDFILE"]

		# ================ GET INSTANCE PARAMETERS AND STORE THEM IN A LIST
		
		fileCount = int(os.environ["PLUGINCONFIG__COUNT"])
		files = []
		for i in range(fileCount):
			item       	= {}
			key        	= "PLUGINCONFIG__ITEMS__"+str(i)
			barcode    	= os.environ[key+"__BARCODE"]
			sample	   	= os.environ[key+"__SAMPLE"]
			input 		= self.analysisDir +"/" + barcode + "_rawlib.bam"
			
			sample = sample.replace(' ', '_')
			
			item["sample"] 	= sample
			item["barcode"] = barcode
			item["input"] 	= input
			
			files.append(item)
			
		# ================ Create output file and structures
		
		htmlOut = open('plotCoverage.html', 'w')
		amplicon_cov_for_sample = {}
		amplicon_cov = {}
		samples_list = []
		
		# ================ LOOP ON EACH FILES AND START COMPUTATION
		
		for item in files:
			print "Analyse {sample}".format(sample=item["sample"])
			
			# ======= utils
			
			previous_chr = None
			sample_name = item["sample"]
			samples_list.append(sample_name)
			amplicon_lst = []
			
			# =======  Get coverage per amplicon
			
			out = os.path.join(self.outputDir, item["sample"]+".cov")
			
			targetReadCoverage = "{pluginDir}/targetReadCoverage.pl -u -a {bam} {bed} | sort +0.3n -1 +1n -2 +2n -3 ".format(pluginDir=self.pluginDir,bam=item["input"], bed=self.bed)
			p = subprocess.Popen(targetReadCoverage, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
			stdout, stderr = p.communicate()
			if p.returncode != 0:
				raise Exception(stderr)
			else:
				for line in stdout.splitlines():
					if line.startswith('chr'):
						amplicon_id = line.split()[3]
						total_reads = line.split()[9]
						total_reads = log10(float(total_reads))
						total_reads = str(total_reads)
						
						#retrieve chromosomes list
						chr_id = line.split()[0]
						if (previous_chr is None):
							previous_chr = chr_id
						else:
							if (chr_id == previous_chr):
								next
							else:
								previous_chr = chr_id
								
						if sample_name in amplicon_cov_for_sample.keys():
							amplicon_cov_for_sample[sample_name] = amplicon_cov_for_sample[sample_name] + "\t" + total_reads
						else:
							amplicon_cov_for_sample[sample_name] = total_reads
						
						if amplicon_id in amplicon_cov.keys():
							amplicon_cov[amplicon_id] = amplicon_cov[amplicon_id] + "\t" + total_reads
						else:
							amplicon_cov[amplicon_id] = total_reads
						
						amplicon_lst.append(amplicon_id)
		
		# print "Amplicon_cov_for_sample\n"
		# for key, value in amplicon_cov_for_sample.items():
		# 	print key, value
		# 	
		# print "Amplicon_cov\n"
		# for key, value in amplicon_cov.items():
		# 	print key, value
		# 	
		# print "Amplicon_lst"
		# print amplicon_lst

		#=== list of files containing amplicon log10(total reads) : one file with all amplicons and 4 others (subset of amplicon to have less data in 1 plot to interpret easier)
		files_list = []
		
		#=== tab file containing header=region_id and lines=samples amplicon total reads for all amplicons and all samples
		all_cov_path = os.path.join(self.outputDir, "all_amplicon_samples.cov.tsv")
		all_cov_file = open(all_cov_path, 'w')
		#print header = amplicon id
		all_cov_file.write("sample\t"+'\t'.join(amplicon_lst)+"\n")
		for sample in samples_list:
			all_cov_file.write(sample + "\t" + amplicon_cov_for_sample[sample] + "\n")
		all_cov_file.close()
		files_list.append(all_cov_path)
		
		#=== convert column to row to make file usable by biologist
		all_infile = open(all_cov_path,'r')
		csv_all_infile = csv.reader(all_infile,delimiter='\t')
		lines = []
		for line in csv_all_infile:
			lines.append(line)
		all_infile.close()
		
		#=== transpose the data and go to log10
		log10_cov_file = os.path.join(self.outputDir, "all_samples_amplicons_log10.cov.tsv")
		log10_out_file = open(log10_cov_file,'w')
		
		transposed = [[lines[j][i] for j in range(len(lines))] for i in range(len(lines[0]))]
		for i in range(len(transposed)):
			for j in range(len(transposed[i])):
				if j==0:
					log10_out_file.write(transposed[i][j])
				else:
					line = "\t"+transposed[i][j]
					log10_out_file.write(line)
			log10_out_file.write("\n")
		log10_out_file.close()
		cmd = "sed -i 's/\./,/g' "+log10_cov_file
		os.system(cmd)
		
		#=== generate amplicon and subset files
		if len(samples_list)>10:
			#subset data files all amplicov cov for 8 samples max
			chunks = [samples_list[i:i+8] for i in xrange(0,len(samples_list),8)]
			sample_subset = 0
			for chunk in chunks:
				sample_subset = sample_subset + 1
				sample_subset_path = os.path.join(self.outputDir, "all_amplicon-sample_subset"+str(sample_subset)+".tsv")
				sample_subset_file = open(sample_subset_path,'w')
				sample_subset_file.write("sample\t" + '\t'.join(amplicon_lst) + "\n")
				for sample in chunk:
					sample_subset_file.write(sample + "\t" + amplicon_cov_for_sample[sample] + "\n")
				sample_subset_file.close()
				files_list.append(sample_subset_path)
				
		#=== 4 subset graphs to view less amplicon info in one graph for all samples
		cmd  = "head -n 1 "+ all_cov_path +"| wc -w"
		nb_columns = os.popen(cmd).read()
		parts = 4
		nb_col_subset = int((int(nb_columns)-1)/parts)
		
		start = 1
		end = start + nb_col_subset - 1
		for i in range(1,5):
			if (end>nb_columns):
				end = nb_columns
			cmd = "cut -f 1,"+str(start)+"-"+str(end)+" "+all_cov_path+" > "+self.outputDir+"/amplicon_subset"+str(i)+"-all_samples.tsv"
			#print cmd + "\n"
			os.system(cmd)
			subset_file = self.outputDir+"/amplicon_subset"+str(i)+"-all_samples.tsv"
			files_list.append(subset_file)
			
			if len(samples_list)>10:
				#subset files containing 8 samples max and amplicons_number/4 (mandatory by IH)
				sample_subset = 0
				for chunk in chunks:
					sample_subset = sample_subset + 1
					sample_subset_path = self.outputDir+"/all_amplicon-sample_subset"+str(sample_subset)+".tsv"
					cmd = "cut -f 1,"+str(start)+"-"+str(end)+" "+sample_subset_path+" > "+self.outputDir+"/amplicon_S"+str(i)+"-samples_S"+str(sample_subset)+".tsv"
					#print cmd + "\n"
					os.system(cmd)
					sub_file_path = self.outputDir+"/amplicon_S"+str(i)+"-samples_S"+str(sample_subset)+".tsv"
					files_list.append(sub_file_path)
					
			start = start + nb_col_subset
			end = end + nb_col_subset
						
		#=== htmlOut

		htmlOut.write('''
		<?xml version="1.0" encoding="iso-8859-1"?>
		<!DOCTYPE HTML>
		<html>
		<!-- <base target="_parent"/> -->
		<META HTTP-EQUIV="Content-Type" CONTENT="text/html; charset=utf-8">
		<head>
		<meta charset="UTF-8">
		<meta name="viewport" content="width=device-width, initial-scale=1.0">
	
		<link rel="stylesheet" media="all" href="/site_media/resources/bootstrap/css/bootstrap.min.css" />
		<link href="/site_media/resources/kendo/styles/kendo.common.min.css" rel="stylesheet" />
		<link href="/site_media/resources/less/kendo.tb.min.css" rel="stylesheet" />
		<link type="text/css" rel="stylesheet" href="/site_media/resources/styles/tb-layout.css" />
		<link type="text/css" rel="stylesheet" href="/site_media/resources/styles/tb-styles.min.css" />
	
		<link rel="stylesheet" type="text/css" href="/site_media/stylesheet.css"/>
		<link rel="stylesheet" type="text/css" href="/site_media/resources/styles/print.css" media="print" />
		<link rel="stylesheet" type="text/css" href="/site_media/resources/styles/report.css" media="screen" />
	
		<script type="text/javascript" src="/site_media/resources/jquery/jquery-1.8.2.min.js"></script>
		<script type="text/javascript" src="/site_media/resources/scripts/kendo.custom.min.js"></script>
	
		<style type="text/css">
		  body {background:white}
		  .help {cursor:help; border-bottom: 1px dotted #A9A9A9}
		</style>
		</head>
	
		<title>Torrent Coverage Plot Report</title>
		<body>
	
		<div class="container-fluid">
	
		<h1><center>Coverage Plot Report</center></h1>
		<h3><center>%s<center></h3>
		<h2><center>Liste des figures</center></h2>
		''' % self.envDict['TSP_ANALYSIS_NAME'])
	
		#=== create amplicons cov graphes
		cpt = 1
		for f in files_list:
			print f + "\n"
			#figure(figsize=(150,50)
			fname = f.split("/")[-1]
			prefix = fname.split(".")[0]
			print "FILE NAME = "+fname+"\n"
			print "PREFIX = "+prefix+"\n"
			
			cov_file = open(f,'r')
			
			my_colors = ['blue','green','tomato','grey','sienna','orange','palevioletred','darkorchid','lime','magenta',
			       'dodgerblue','lightcoral','lightblue','darkkhaki','cyan','yellow','plum','peru','steelblue',
			       'mediumspringgreen','red','yellowgreen','maroon','gold','black','darkseagreen','burlywood',
			       'deeppink','slategrey','greenyellow','indigo','hotpink','firebrick','indianred','mistyrose',
			       'darkolivegreen','olive','pink','orangered','navajowhite','palegreen','darkslategrey','seashell',
			       'fuchsia','papayawhip','blanchedalmond','chartreuse','dimgray','peachpuff','springgreen','aquamarine',
			       'white','lightsalmon','darkslategray','brown','ivory','darkgrey','lawngreen','chocolate','crimson',
			       'forestgreen','slateblue','lightseagreen','mintcream','silver','antiquewhite','mediumorchid',
			       'skyblue','gray','darkturquoise','goldenrod','darkgreen','floralwhite','darkviolet','darkgray',
			       'moccasin','saddlebrown','darkslateblue','lightskyblue','lightpink','mediumvioletred','limegreen',
			       'darkmagenta','palegoldenrod','turquoise','lightgrey','lightgoldenrodyellow','darkgoldenrod','lavender',
			       'sandybrown','thistle','violet','navy','dimgrey','tan','rosybrown','olivedrab','ghostwhite','honeydew',
			       'cornflowerblue','linen','darkblue','powderblue','seagreen','snow','mediumblue','royalblue','lightcyan',
			       'mediumpurple','midnightblue','cornsilk','bisque','slategray','darkcyan','khaki','wheat','teal',
			       'deepskyblue','salmon','darkred','lightslategray','aliceblue','lightslategrey','lightgreen','orchid',
			       'gainsboro','mediumseagreen','mediumturquoise','lemonchiffon','cadetblue','lightyellow','lavenderblush',
			       'coral','purple','aqua','whitesmoke','mediumslateblue','darkorange','mediumaquamarine','darksalmon',
			       'beige','blueviolet','azure','lightsteelblue','oldlace']
			nb = 0
			for line in cov_file:
				if line.startswith('sample'):
					x_names = line.split('\t')
					x_names.pop(0) #remove first element of the list (header line, remove "sample" label and retreive amplicon id for x_axis)
					x_len = len(x_names)
					x_data = range(len(x_names))
					x = array(x_data)
					
				else:
					y_data = line.split('\t')
					y_legend = y_data[0] #first field of the line = sample name for legend
					y_data.pop(0)
					y2_data =[]
					#fig=figure(1,figsize=(50,20))
					fig=figure(1,figsize=(65,20))
					for n in y_data:
						y2_data.append(n)
						
					y2 = array(y2_data)
					my_color = my_colors[nb]
					plot(x,y2,label=y_legend,color=my_color,linewidth=3.5)
					nb = nb + 1
			cov_file.close()
			width = 200
			height = 5000
			for i in range(len(x_names)):
				amp = x_names[i]
			
			xticks(x_data,x_names,rotation=90,fontsize=18)
			yticks(fontsize=14)
			plt.axis('tight')
			ylim(ymin=0,ymax=5)
			xlabel("amplicons",fontsize=20)
			ylabel("total reads (log10)",fontsize=20)
			axhline(y=3,color='black',linestyle='--',label='1000')
			axhline(y=2.7,color='b',linestyle='--',label='500')
			axhline(y=2.48,color='orange',linestyle='--',label='300')
			axhline(y=2,color='r',linestyle='--',label='100')
			#lgd = legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=4, mode="expand", borderaxespad=0.)
			legend(loc='center left',bbox_to_anchor=(1,0.5)) #legend(fontsize=20,loc='center left',bbox_to_anchor=(1,0.5))
			#autoscale(enable=True)
			#tight_layout(rect=[0,0,0.9,1])
			
			#show()
			
			#output file graph
			file_name = prefix+'_plot.png'
			img_path = os.path.join(self.outputDir,file_name)
			# fig.autofmt_xdate() # inclinaison des legendes
			fig.savefig(img_path,dpi=80)
			#fig.savefig(img_path,dpi=100,bbox_inches='tight',bbox_extra_artists=[legend])
			#savefig(img_path)
			htmlOut.write('<h4><b>Figure %s:</b> %s <br/></h4>' % (cpt,file_name.replace(".png","")))
			htmlOut.write('<img src="%s" /> ' % file_name)
			cpt = cpt+1
			
			close()
				
		return True
	
if __name__ == "__main__":
  PluginCLI()
