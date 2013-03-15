#!/usr/bin/python

import csv
import os.path
from collections import namedtuple
import sn
import os
import sys,string
import numpy as np
import math
import vcf
import fnmatch


#try:	
#	file_map = sys.argv[1];dir_files_phenotype1 = sys.argv[2];dir_files_phenotype2 = sys.argv[3];outfilename = sys.argv[4]

#except:
#     print "Usage:",sys.argv[0], "file.map dir_files_phenotype1 dir_files_phenotype2 outfile";sys.exit(1)



file_map="/home/cristovao/Desktop/AUS_project/public_datas/opensnp_datadump.201303070733/phenotypes_201303070733.csv"

folder="/home/cristovao/Desktop/AUS_project/public_datas/opensnp_datadump.201303070733"


def get_dataset():
    """return """
    handle = csv.DictReader(open(file_map, "r"),
            #fieldnames=["user_id","date_of_birth","chrom_sex","Jewish Ancestry","Subjective dream intensity","Webbed toes","Dyslexia","Artistic ability","lips size","ethnicity","Acrophobia","Myers-Briggs Type Indicator","Irritable Bowel Syndrome","Diego Blood Group","Cholesterol","Moles raised","Autism","Interest in Spirituality and Mysticism","Physician-diagnosed celiac/coeliac disease","Hypertriglyceridemia","SAT Writing","Panic Disorder","Bone Mineral Density","Sexual Preferences","Energy Level","Faktor 5 Leiden (F5)","Age learned to read","ear proximity to head ","Atheism","Earwax type","ring finger longer than index finger","Eye with Blue Halo ","Beard Color","Birth year","Migraine frequency","Serotonin transporter","Sport interest","Number of toes","Number of wisdom teeth","Widow's Peak","natural skinny","Wake up preference","Lisp","Do you like the taste of hops?","Wanting to be immortal","Purposefulness ","Ambition","Do hops taste like soap?","ABH Blood Group (antigens) ","Fish Preference","Smell of coffee in urine","hair on fingers","Neanderthal","Are You The Advertising Phenotype?","(male) penis releases pre-cum when sexually aroused.","Morton's Toe","Sports interest","Does cilantro taste like soap to you?","Tongue roller","Enjoy watching TV","Aspirin Allergy","libido ","Blood type","First word","Enjoy using the Internet","mtDNA Haplogroup (PhyloTree)","Like the taste of Stevia","Negative reaction to fluoroquinolone antibiotics","white skin","Fat-pad knee syndrome","Ability to Tan","Strabismus","Amblyopia","Autoimmune disorder","Y-DNA Haplogroup (ISOGG)","Asthma","Freckling","form of the nose","Ancestry","Metabolic Syndrome [MetS]","Enjoy riding a motorbike","Hair Color","Tea consumption","Height","Sex","Motion sickness","Cystic Fibrosis Like Disease","mouth size","Peanut butter preference","Sneezing induced by sexual ideation or orgasm?","Woolnerian Tip (Darwin's Tubercle)","SAT Math","prognathism","Taste of broccoli","Jogger","Phobia","Kell Blood Group (K/k antigens)  ","Desmoid Tumor","SAT Verbal","Astigmatism","excessive daytime sleepiness","Enjoy driving a car","ABO Rh ","Kidd Blood Group","Sense of smell","apthous in mouth tendency","Allergic/bad reaction to fish oil supplements","Interested in news from real newspaper / news from the Internet","erectil disfunction ","Index Toe Longer than Big Toe","Hair Type","Penis Circumference at Glans","Penis Length","Intolerance: gluten, casein, soy","Weight","Short-sightedness (Myopia)","brown hair colour","SAT - when taken","Anorgasmia","Nicotine dependence","CMV serostatus","Musical Perfect Pitch","Rheumatoid Arthritis","(Male) Nipple's size","ADHD","Insect bites and stings","Colour Blindness","Lactose intolerance","Have ME/CFS","Atypical Sulfonomide Antibiotic Reaction","Cramps","Political Ideology","Handedness","cluster headache","Eye color","Social Level","Earlobe:  Free or attached","Photic Sneeze Reflex (Photoptarmis)","Coffee consumption","Penicillin reaction","Do you have a parent who was diagnosed with Alzheimer's disease?","R1b1a2a1a1b","Good / poor eater as child","Abnormal Blood Pressure","Type II Diabetes","Migraine","Colon cancer ONLY FOR (rs3219489 GG)!","Ability to find a bug in openSNP","Eurogenes","head form","Cleverness","ENTP","Can you smell cut-grass?","Asparagus Metabolite Detection"],
        delimiter=";")

    return handle


def get_user(pheno, variation):
    """Return list of the user with a specific variation """
    dataset = get_dataset()
    user_list = []
    for i in dataset:
        if i[pheno] == variation:
            user_list.append(i["user_id"])
    dataset=[]
    return user_list





def create_dir(user_list,variation):
    """Create a folder from a list of the user"""
    user_list=list(set(user_list))
    print "total of the user", len(user_list), user_list
    files= os.listdir(folder)
    
    #variation="_".join(variation.split())
    
    os.system("mkdir "+variation)
    n=0
    for j in user_list: 
        for i in files:
          if fnmatch.fnmatch(i, '*.txt'):   
            u="user"+j+"_"
            if u in i:
                print i
                os.system("cp "+folder+"/"+i +" " +variation+"/")
                n=1+n
    print "total of the files copied", n        


#------------------  execution  ------------------- "Eye color"

fieldnames=open(file_map).readline().split(';')
fieldnames.sort()
print "\n\n---------------------------  fieldnames (Phenotypes)\n"


for i in fieldnames:
    print i



p=raw_input("\n---------------------------  Phenotype: ")
variations_list=[] 
for i in get_dataset():
    if not i[p] in variations_list: 
        variations_list.append(i[p])
        print i[p]        


v=raw_input("\n---------------------------  Variations: ")
v=v.split(";")

print "\n"

os.system("mkdir "+"_".join(p.split()))


for i in v:
    print "Variations: ", i
    l=get_user( p, i)
    variation="_".join(i.split())
    create_dir(l,variation)
    os.system("mv "+ variation+" "+"_".join(p.split()))    
    print "\n"











