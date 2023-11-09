import math
import random
import itertools
import operator

ln2 = math.log(2)

##############
#Parameters. #
##############
#immune shape space parameters
gene_len = 20  # see Smith and Perelson PNAS 1999
gene_vocab = 4  # see Smith and Perelson PNAS 1999
max_dist = 7
min_dist = 4
germline_affinity = 6
isotype_position = gene_len + 1
#immune system parameters
differentiation_rate = 0.10  # see (1 - recycling rate) Oprea & Perelson J. Immunology 1997
mutation_rate = 0.10  # see Oprea & Perelson J. Immunology 1997
isotype_rate = 0.0
lethal_fraction = 0.0

######################################
#ANTIGEN BINDING AFFINITY FUNCTIONS. #
######################################
def BindingAffinity(phenotype, affinity_factor):
    #right now binding affinity scales by ~1 with frequency of phenotype
    if (phenotype <= min_dist):
        return math.pow(affinity_factor, germline_affinity - (float(min_dist)))
    elif (phenotype <= max_dist):
        return math.pow(affinity_factor, germline_affinity - (float(phenotype)))
    else:
        return 0.0
class BindingAffinityPrecalc:
    def __init__(self):
        self._BA_list = []
        self._aff = []

    def add_aff(self, value):
        BA_array = []
        for i in range(0, max_dist + 1):
            BA_array.append(BindingAffinity(i, value))
        self._aff.append(value)
        self._BA_list.append(BA_array)
        print("ADDING AFF VALUE " + str(value))

    def aff_value(self, value):
        try:
            i = self._aff.index(value)
            return i
        except ValueError:
            self.add_aff(value)
            return self._aff.index(value)

    def get_BA(self, phenotype, affinity_factor):
        if (phenotype <= max_dist):
            BA_array = self._BA_list[self.aff_value(affinity_factor)]
            return BA_array[phenotype]
        else:
            return 0.0
BA = BindingAffinityPrecalc()

def Crossreactivity( ep1, ep2 ):
	ne = operator.ne
	hamming_dist = sum(map(ne, ep1.get_sequence(), ep2.get_sequence()))
	cross_matrix = [1.0, 0.35, 0.27, 0.16, 0.11, 0.070, 0.044, 0.027]
	cross = 0
	if (hamming_dist < len(cross_matrix)):
		cross = cross_matrix[ hamming_dist ]
	return cross 

def AvgCrossreactivity( antigen_list ):
	n_antigen = len(antigen_list)
	n_epitope = len(antigen_list[0].epitope_all())
	cross = []
	for i in range(0, n_epitope):
		cross.append(float(0))
		for j in range(0, n_antigen):
			for k in range(0, n_antigen):
				if (j != k ):
					ep1 = antigen_list[j].epitope(i)
					ep2 = antigen_list[k].epitope(i)
					cross[i] = cross[i] + Crossreactivity(ep1, ep2)/(n_antigen-1)
	for i in range(0, len(cross)):
		cross[i] = cross[i] / float(n_antigen)
	avg = sum(cross)/float(n_epitope)
	return avg

def BindingAffinity_Pre(phenotype, affinity_factor):
    return BA.get_BA(phenotype, affinity_factor)

def ApparentSize(value, phenotype, aff_factor):
    if (float(phenotype) <= max_dist):
        n = float(value) * BindingAffinity_Pre(phenotype, aff_factor)
        return n
    else:
        return 0.0

######################
#BCR GENE FUNCTIONS. #
######################
def GeneRandom():
    sequence = ""
    for i in range(0, gene_len):
        sequence += str(random.randint(1, gene_vocab))
    sequence += 'M'  #add isotype
    return sequence

def GeneMutate(sequence):
    mutation_position = random.randint(1, gene_len)
    out_seq = ""
    for i in range(0, gene_len):
        if ( i == mutation_position ):
            out_seq += str(random.randint(1, gene_vocab))
        else:
            out_seq += sequence[i]
    out_seq += sequence[isotype_position - 1]  #add isotype
    return out_seq

def IsotypeSwitch(sequence):
    out_seq = ""
    for i in range(0, gene_len):
        out_seq += sequence[i]
    out_seq += 'G'  #add isotype
    return out_seq

def GenePhenotype(sequence, antigen_in):
    ne = operator.ne
    epitope = antigen_in.epitope_all()
    dist = 999
    for i in range(0, len(epitope)):
        hamming_dist = sum(map(ne, sequence, epitope[i].get_sequence()))
        if (hamming_dist < dist ):
            dist = hamming_dist
    phenotype = min(max_dist + 1, dist)
    return phenotype

def GeneEpitope(sequence, antigen_list):
    ne = operator.ne
    dist = 999
    epitope_num = 999
    for antigen in antigen_list:
        epitope = antigen.epitope_all()
        for i in range(0, len(epitope)):
            hamming_dist = sum(map(ne, sequence, epitope[i].get_sequence()))
            if (hamming_dist < dist ):
                dist = hamming_dist
                epitope_num = i
    return epitope_num

def GeneFromPhenotype(phenotype, antigen, ep_num):
    sequence = ""
    rand_pos = []
    for i in range(1, gene_len + 1):
        rand_pos.append(i)
    random.shuffle(rand_pos)
    num_match = gene_len - phenotype
    epitope = antigen.epitope(ep_num).get_sequence()
    for i in range(0, gene_len):
        if (rand_pos[i] <= num_match):
            sequence += epitope[i]
        else:
            new_val = epitope[i]
            while ( new_val == epitope[i] ):
                new_val = str(random.randint(1, gene_vocab))
            sequence += new_val
    sequence += 'M'  #add isotype
    return sequence

#################################
#IMMUNE SYSTEM COMPONENT TYPES. #
#################################
#Name: Population
#Desc: A base class that describes the  population of a given immune system component.
#It is made up of subpopulations, each with their own genotype and size. There are two types of Populations, Antigen and BCell (see below).
class Population:
    def __init__(self, name, value):
        self._name = name
        self._value = value
        self._type = ''

    def set_type(self, type):
        self._type = type

    def set_size(self, value):
        self._value = value

    def size(self):
        return self._value

    def decrease(self, n):
        self._value = self._value - n

    def increase(self, n):
        self._value = self._value + n

    def name(self):
        return str(self._name)

#Name: GroupPopulation
#Desc: An aggregate population that is made up of individual Population objects
class GroupPopulation:
    def __init__(self, name):
        self._name = name
        self._pop = []
        self._value = 0

    def add_population(self, pop):
        self._pop.append(pop)

    def calc_population(self):
        self._value = 0
        for i in range(0, len(self._pop)):
            self._value += self._pop[i].size()

    def size(self):
        self.calc_population()
        return self._value

    def length(self):
        return (len(self._pop))

    def random_select(self):
        self.calc_population()
        r = random.random()
        select = 0.0
        for i in range(0, len(self._pop)):
            select += float(self._pop[i].size()) / float(self._value)
            if (r <= select ):
                return i
                break

    def decrease(self, n):
        self._pop[self.random_select()].decrease(n)

    def sub_decrease(self, i, n):
        self._pop[i].decrease(n)

    def increase(self, n):
        self._pop[self.random_select()].increase(n)

    def sub_increase(self, i, n):
        self._pop[i].increase(n)

#Name: Epitope
#Desc: Defines a single epitope, which includes an epitope sequence (in immune shape space), as well as its relative immunogenicity and clearance rates
class Epitope:
    def __init__(self, name, sequence, immunogenicity, clearance):
        self._name = name
        self._sequence = sequence
        self._immunogenicity = immunogenicity
        self._clearance = clearance

    def get_sequence(self):
        return self._sequence

    def immunogenicity(self):
        return self._immunogenicity

    def clearance(self):
        return self._clearance

#Name: Antigen
#Desc: Defines an antigen population, which is one type of immune system component. It is defined by a name, a population size, and an AntigenType
class Antigen(Population):
    def __init__(self, name, value, antigen):
        """ :rtype : object """
        Population.__init__(self, name, value)
        self._antigen = antigen

    def change_antigen(self, new_antigen):
        self._antigen = new_antigen

    def return_antigen(self):
        return self._antigen

#Name: AntigenType
#Desc: Defines an antigen type (such as a given antigen strain). This includes a name, a numerical antigenID, and a list of epitopes that make up the antigen
#note: in future version, antigenID should be set internally for ease of use
class AntigenType:
    def __init__(self, name, antigenID):
        self._name = name
        self._epitope_array = []
        self._antigenID = antigenID

    def add_epitope(self, ep):
        self._epitope_array.append(ep)

    def epitope(self, r):
        return self._epitope_array[r]

    def epitope_all(self):
        return self._epitope_array

    def epitope_num(self):
        return len(self._epitope_array)

    def reset_epitopes(self):
        self._epitope_array = []

    def ID(self):
        return self._antigenID

#Name: Bcell 
#Desc: Defines a B cell population. It is indexed by genotype, and keeps track of the genotype, the size of the genotype population, 
# the Epitope that genotype recognizes, and the phenotype of that genotype with respect to all Antigens in the system.
class BCell( Population ):
	def add_genotype( self, gene, n ):
		self._genotype_list.append( gene )
		epitope = GeneEpitope(gene, self._antigen_list)
		self._epitope_list.append( epitope )
		self._size_list.append( n )
		for antigen in self._antigen_list:
			phenotype = GenePhenotype(gene, antigen)
			phenotype_list = self._phenotype_master[antigen.ID()]
			phenotype_list.append( phenotype )
			subpopulation = self._subpopulation_master[antigen.ID()]
			subpopulation[epitope][phenotype] += n

	def phenotype_size( self, phenotype, antigenID ):
		pop_size = 0
		phenotype_list = self._phenotype_master[antigenID]
		for i in range(0, len(phenotype_list)):
			if (phenotype_list[i] == phenotype):
				pop_size += self._size_list[i]
		return pop_size

	def epitope_size( self, epitope ):
		pop_size = 0
		for i in range(0, len(self._epitope_list)):
			if (self._epitope_list[i] == epitope):
				pop_size += self._size_list[i]
		return pop_size

	def RealSizePhenotype( self, epitope, phenotype, type, antigenID):
		subpopulation = self._subpopulation_master[antigenID]
		n = subpopulation[epitope][phenotype]
		alpha = 0.0
		if (type == "immunogenicity"):
			alpha = self._antigen_list[antigenID].epitope( epitope ).immunogenicity()
		elif (type == "clearance"):
			alpha = self._antigen_list[antigenID].epitope( epitope ).clearance()
		else:
			print("ERROR in APPARENTSIZEPHENOTYPE")
		real_size = n * alpha
		return real_size

	def ApparentSizePhenotype( self, epitope, phenotype, aff_factor, type, antigenID):
		subpopulation = self._subpopulation_master[antigenID]
		n = subpopulation[epitope][phenotype]
		alpha = 0.0
		if (type == "immunogenicity"):
			alpha = self._antigen_list[antigenID].epitope( epitope ).immunogenicity()
		elif (type == "clearance"):
			alpha = self._antigen_list[antigenID].epitope( epitope ).clearance()
		else:
			print("ERROR in APPARENTSIZEPHENOTYPE")
		app_size = ApparentSize( n, phenotype, aff_factor ) * alpha
		return app_size

	def ApparentSizeEpitope( self, epitope, aff_factor, type, antigenID):
		app_size = 0.0
		for i in range(0,8):
			app_size += self.ApparentSizePhenotype( epitope, i, aff_factor, type, antigenID )
		return app_size

	def ApparentSizeAll( self, aff_factor, type, virus ):
		antigenID = virus.return_antigen().ID()
		app_size = 0.0
		for i in range(0, self._antigen_list[0].epitope_num() ):
			app_size += self.ApparentSizeEpitope( i, aff_factor, type, antigenID)
		return app_size

	def select_random( self ):
		r = random.random()
		outcome = float(0.0)
		for i in range(0, len(self._genotype_list)):
			outcome += float(self._size_list[i])/float(Population.size( self ))
			if (r <= outcome):
				return str(self._genotype_list[i])
				break
		print(("ERROR: COULD NOT FIND RANDOM select random " + str(self.name()) ))

	def select_epitope_phenotype( self, epitope, phenotype, antigenID ):
		r = random.random()
		outcome = float(0.0)
		subpopulation = self._subpopulation_master[antigenID]
		phenotype_list = self._phenotype_master[antigenID]
		pop_size = subpopulation[epitope][phenotype]
		tot_size = 0
		for i in range(0, len(self._genotype_list)):
			#correct = 0.0
			if (phenotype_list[i] == phenotype and self._epitope_list[i] == epitope ):
				#correct = 1.0
				tot_size += self._size_list[i]
				outcome += (float(self._size_list[i]) )/float(pop_size)
			if (r <= outcome):
				return str(self._genotype_list[i])
				break
		print(("ERROR: COULD NOT FIND ###RANDOM select phenotype " + str(self._name) + ' ' + str(r) + ' ' +str(outcome) + ' ' + ' ' + str(tot_size) + ' ' + str(pop_size)))

	def select_random_weighted( self, aff_factor, max_rate, type, agg_rate, antigen ):
		r = random.random()
		base_phenotype = 999
		base_epitope = 999
		total_p = 0.0
		for i in range(0, antigen.epitope_num() ):
			for j in range(0,8):
				p_size = min(self.ApparentSizePhenotype( i, j, aff_factor, type, antigen.ID())*agg_rate , self.RealSizePhenotype(i, j, type, antigen.ID())*max_rate)
				total_p += p_size
		found_it = 1
		outcome = 0.0
		for i in range(0, antigen.epitope_num() ):
			for j in range(0,8):
				p_size = min(self.ApparentSizePhenotype( i, j, aff_factor, type, antigen.ID())*agg_rate , self.RealSizePhenotype(i, j, type, antigen.ID())*max_rate)
				outcome += p_size/total_p
				if (r <= outcome  and found_it == 1):
					base_epitope = i
					base_phenotype = j
					found_it = 0

		if ( base_epitope == 999 or base_phenotype == 999):
			print("ERROR PHENOTYPE/EPITOPE NOT SELECTED IN REACT")
		subpopulation = self._subpopulation_master[antigen.ID()]
		pop_size = subpopulation[base_epitope][base_phenotype]
		if (pop_size <= 0):
			print(("ERROR: POPULATION SIZE FOR #PHENOTYPE " + str(base_phenotype) + " IS ZERO " + self.name() ))
		return self.select_epitope_phenotype( base_epitope, base_phenotype, antigen.ID() )

	def diversity(self, threshold):
		line_count = 0
		for i in range(0, len(self._genotype_list)):
			if (self._size_list[i] >= threshold):
				line_count = line_count + 1
		return line_count

	def diversity2( self, threshold):
		num = threshold * float(Population.size( self ))
		sort_list = sorted( self._size_list, reverse = True )
		gene_count = 0
		pop_count = 0
		for i in range(0, len(self._genotype_list)):
			if (pop_count < num):
				pop_count += sort_list[i]
				gene_count += 1
		return gene_count

	def genotype_increase( self, gene, n):
		Population.increase( self, n )
		try:
			i = self._genotype_list.index( gene )
			self._size_list[i] += n
			for antigen in self._antigen_list:
				subpopulation = self._subpopulation_master[antigen.ID()]
				phenotype_list = self._phenotype_master[antigen.ID()]
				subpopulation[self._epitope_list[i]][phenotype_list[i]] += n
		except ValueError:
			self.add_genotype( gene, n )

	def genotype_decrease( self, gene, n ):
		Population.decrease(self, n )
		try:
			i = self._genotype_list.index( gene )
			for antigen in self._antigen_list:
				phenotype_list = self._phenotype_master[antigen.ID()]
				subpopulation = self._subpopulation_master[antigen.ID()]
				subpopulation[self._epitope_list[i]][phenotype_list[i]] -= n
			if ((self._size_list[i] - n) <= 0):
				self._size_list.pop(i)
				self._genotype_list.pop(i)
				self._epitope_list.pop(i)
				for antigen in self._antigen_list:
					phenotype_list = self._phenotype_master[antigen.ID()]
					phenotype_list.pop(i)
			else:
				self._size_list[i] -= n
				
		except ValueError:
			print((str(self._name) + " "+ str(Population.size(self)) +str(gene)+" is not in the list!!"))

	def increase( self, n ):
		phenotype = 7
		r = random.random()
		#if (r <= 0.10):
		#	phenotype = 6
		#if (r <= 0.01):
		#	phenotype = 5
		n_antigen = len(self._antigen_list)
		n_epitope = len(self._antigen_list[0].epitope_all())
		r_antigen = random.random()
		r_epitope = random.random()
		count = 0.0
		antigen = 999
		for i in range(0,n_antigen):
			count = count + float(1.0/n_antigen)
			if (r_antigen <= count and antigen == 999):
				antigen = i
		epitope = 999
###NEW###
#		if (r_epitope <= float(1.0/(n_antigen+1))):
#			epitope = 0
#		else:
#			epitope = 1
		cross = []
		for i in range(0, n_epitope):
			cross.append(float(0))
			for j in range(0, n_antigen):	
				if (j != antigen ):
					ep1 = self._antigen_list[antigen].epitope(i)
					ep2 = self._antigen_list[j].epitope(i)
					cross[i] = cross[i] + Crossreactivity(ep1, ep2)/(n_antigen-1)
		correction_factor = []
		for i in range(0,n_epitope):
			correction_factor.append( float(1/(1+(n_antigen-1)*cross[i])))
		avg_correction_factor = sum(correction_factor)/len(correction_factor)
		r_epitope = r_epitope * avg_correction_factor
		count = 0.0
		for i in range(0,n_epitope):
			count = count + (float(1.0/n_epitope) * (correction_factor[i]))
			if (r_epitope <= count and epitope == 999):
				epitope = i

		base_gene = GeneFromPhenotype( phenotype, self._antigen_list[antigen], epitope )
		self.genotype_increase( base_gene, n )

	def decrease( self, n ):
		self.genotype_decrease( self.select_random(), 1 )

	def calc_crossreactivity( self, antigen1, antigen2 ):
		cross_num = 0
		spec1_num = 0
		spec2_num = 0
		for i in range(0, len(self._genotype_list)):
			phenotype1 = GenePhenotype( self._genotype_list[i], antigen1 )
			phenotype2 = GenePhenotype( self._genotype_list[i], antigen2 )
			if (phenotype1 <= max_dist and phenotype2 <= max_dist):
				cross_num += self._size_list[i]
			elif (phenotype1 <= max_dist):
				spec1_num += self._size_list[i]
			elif (phenotype2 <= max_dist):
				spec2_num +=self._size_list[i]
		output = [cross_num, spec1_num, spec2_num]
		return output

	def calc_transcend( self, antigen_list ):
		output = []
		for i in range(0,len(antigen_list)+1):
			output.append(0)
		for i in range(0, len(self._genotype_list)):
			strain_num = 0
			#for antigen in self._antigen_list:
			for j in range(0, len(antigen_list)-1):
				phen = GenePhenotype( self._genotype_list[i], antigen_list[j] )
				if (phen <= max_dist):
					strain_num += 1
			output[strain_num] += self._size_list[i]
		return output

	def calc_neutralization( self, antigen1, antigen2 ):
		spec1_num = 0
		spec2_num = 0
		for i in range(0, len(self._genotype_list)):
			phenotype1 = GenePhenotype( self._genotype_list[i], antigen1 )
			phenotype2 = GenePhenotype( self._genotype_list[i], antigen2 )
			if (phenotype1 <= max_dist):
				spec1_num += self._size_list[i]*BindingAffinity_Pre(phenotype1, 2.5)
			if (phenotype2 <= max_dist):
				spec2_num +=self._size_list[i]*BindingAffinity_Pre(phenotype2, 2.5)
		output = [spec1_num, spec2_num]
		return output

	def calc_depletion( self, antigen1, antigen2 ):
		#deplete with antigen2, then measure activity against antigen1
		spec1_num = 0
		spec2_num = 0
		for i in range(0, len(self._genotype_list)):
			phenotype1 = GenePhenotype( self._genotype_list[i], antigen1 )
			phenotype2 = GenePhenotype( self._genotype_list[i], antigen2 )
			if (phenotype1 <= max_dist and phenotype2 > max_dist):
				spec1_num += self._size_list[i]*BindingAffinity_Pre(phenotype1, 2.5)
		output = spec1_num
		return output

	def calc_cross(self, antigen1, epitope):
		subpopulation = self._subpopulation_master[antigen1.ID()]
		pop = [0,0,0,0]
		for phenotype in range(0, 5):
			pop[0] = pop[0] + subpopulation[epitope][phenotype]
		pop[1] = subpopulation[epitope][5]
		pop[2] = subpopulation[epitope][6]
		pop[3] = subpopulation[epitope][7]
		return pop

	def calc_isotype( self ):
		g_num = 0
		for i in range(0, len(self._genotype_list)):
			gene = self._genotype_list[i]
			if (gene[isotype_position - 1] == 'G'):
				g_num += self._size_list[i]
		m_num = self.size() - g_num
		output = [m_num, g_num]
		return output

	def generate_population( self, value):
		for antigen in self._antigen_list:
			n_epitope = len(antigen.epitope_all())
			for j in range(0, n_epitope):
				for phenotype in range(7,8):
					current_pop = 0
					pop = int(math.floor(math.pow(10,-1*(5+max_dist-phenotype)) * float(value))/n_epitope)
					while (current_pop < pop):
						current_pop += 1
						base_gene = GeneFromPhenotype( phenotype, antigen, j )
						min_phen = 7
						for k in range(0,len(self._antigen_list)):
							phenotype1 = GenePhenotype( base_gene, self._antigen_list[k] )
							if (phenotype1 < min_phen ):
								min_phen = phenotype1
						if (min_phen == 7):
							self.genotype_increase(base_gene, 1)
						print( base_gene )
		high_pop = int(math.floor(math.pow(10,-1*(5+max_dist-phenotype)) * float(value))) * len(self._antigen_list) * 2
		for i in range(0, high_pop):
			for antigen in self._antigen_list:
				n_epitope = len(antigen.epitope_all())
				for j in range(0, n_epitope):
					for phenotype in range(7,8):
						pop = int(math.floor(math.pow(10,-1*(5+max_dist-phenotype)) * float(value))/n_epitope)
						subpopulation = self._subpopulation_master[ antigen.ID() ]
						current_pop = subpopulation[j][phenotype]
						if ( subpopulation[j][phenotype] > pop ):
							del_gene = self.select_epitope_phenotype(j, phenotype, antigen.ID())
							self.genotype_decrease( del_gene, 1)

	def return_gene_list( self ):
		return self._genotype_list

	def return_size_list( self ):
		return self._size_list

	def __init__( self, name, value, antigen_list ):
		Population.__init__( self, name, 0 )
		Population.set_type(self, "bcell")
		self._genotype_list = []
		self._epitope_list = []
		self._size_list = []
		self._antigen_list = antigen_list
		epitope_num = antigen_list[0].epitope_num()
		self._phenotype_master = []
		self._subpopulation_master = []
		for antigen in antigen_list:
			phenotype_list = []
			subpopulation = [[0]*20 for x in range(epitope_num)]#change from epitope_num to 6
			self._phenotype_master.append(phenotype_list)
			self._subpopulation_master.append(subpopulation)
		self.generate_population( value )

###################
#REACTIONS TYPES. #
###################
#OrderZero: Defines the base class for zero-order reactions 
class OrderZero:
	def __init__( self, name, k ):
		self._name = name
		self._k = ln2 * float( k )
		self.__rate = 0

	def rate( self ):
		self.__rate = self._k
		return self.__rate

#OrderOne: Defines the base class for first-order reactions
class OrderOne:
	def __init__(self, name, k, A):
		self._name = name
		self._k = ln2 * float(k)
		self._A = A

	def rate(self):
		self.__rate = self._k * self._A.size()
		return self.__rate

#OrderTwoPhenotype: Defines the base class for second-order reactions where the reaction rate is the result of a heterogenious population, with varying phenotypes (such as B cell stimulation)
class OrderTwoPhenotype:
    def __init__(self, name, k, A, B, max_rate, aff_factor, type):
        self._name = name
        self._k = ln2 * float(k)
        self._A = A
        self._B = B
        self._max_rate = ln2 * max_rate
        self._aff_factor = aff_factor
        BA.aff_value(aff_factor)
        self._type = type

    def rate(self):
        self.__rate = 0.0
        if (self._B.size() > 0):
            Bcell_app = self._B.ApparentSizeAll(self._aff_factor, self._type, self._A)
            self.__rate = min(self._A.size() * Bcell_app * self._k, self._B.size() * self._max_rate)
        return self.__rate

###########################
#IMMUNE SYSTEM REACTIONS. #
###########################
#Stimulation: Defines the reaction for B cell stimulation
class Stimulation(OrderTwoPhenotype):
    def __init__(self, name, k, A, B, C, max_rate, aff_factor, type):
        OrderTwoPhenotype.__init__(self, name, k, A, B, max_rate, aff_factor, type)
        self.A = A
        self.B = B
        self.C = C
        self.aff_factor = aff_factor
        self.type = type
        self.max_rate = max_rate
        self.k = k

    def react(self):
        antigen = self.A.return_antigen()
        agg_rate = self.k * self.A.size()
        base_gene = self.B.select_random_weighted(self.aff_factor, self.max_rate, self.type, agg_rate, antigen)
        self.B.genotype_decrease(base_gene, 1)
        self.C.genotype_increase(base_gene, 1)

#Formation: Defines the spontaneous formation of a B cell
class Formation(OrderZero):
    def __init__(self, name, k, A):
        OrderZero.__init__(self, name, k)
        self.A = A

    def react(self):
        self.A.increase(1)

#Decay: Define the first order decay of many components, such as antigen, antibodies, plasma cells
class Decay(OrderOne):
    def __init__(self, name, k, A):
        self.A = A
        OrderOne.__init__(self, name, k, A)

    def react(self):
        self.A.decrease(1)

#PopulationDecay: Defines the first order decay of a population under a carrying capacity, for example GC B cells
class PopulationDecay:
    def __init__(self, name, r, k, capacity, A, B):
        self.name = name
        self.r = r  #replication rate
        self.k = k  #minimum decay rate
        self.A = A
        self.B = B
        self.capacity = float(capacity)  #carrying capacity
        self.__rate = 0

    def rate(self):
        self.__rate = 0
        population_size = float(self.A.size())
        k_adjust = self.r * (population_size / self.capacity)
        k = max(self.k, k_adjust)
        self.__rate = k * population_size
        return self.__rate

    def react(self):
        #self.A.decrease(1)
        if (self.B.size() > 1):
            self.B.decrease(1)

#Differentiation: Defines the first order reaction of differentiation into one of three components, for example, a stimulated B cell can divide into a daughter B cell, memory cell, or plasma cell.
#The division into a daughter B cell has a probability including a somatic mutation.
class Differentiation(OrderOne):
    def __init__(self, name, k, A, B, C, D, E, max_rate):
        OrderOne.__init__(self, name, k, A)
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.E = E
        self.max_rate = max_rate

    def react(self):
        gene_A = self.A.select_random()
        self.A.genotype_decrease(gene_A, 1)
        iso_rate = isotype_rate
        mut_rate = (mutation_rate + mutation_rate / gene_vocab) * (1 - lethal_fraction)
        diff_rate = differentiation_rate
        for i in range(2):
            new_gene = gene_A
            r = random.random()
            if ( r <= iso_rate ):
                new_gene = IsotypeSwitch(gene_A)
                self.B.genotype_increase(new_gene, 1)
            elif ( r <= mut_rate + iso_rate ):
                new_gene = GeneMutate(new_gene)
                self.B.genotype_increase(new_gene, 1)
            elif (r <= mut_rate + iso_rate + diff_rate):
                if (random.random() < 0.50):
                    self.C.genotype_increase(gene_A, 1)
                elif (random.random() < 0.75):
                    self.D.genotype_increase(gene_A, 1)
                else:
                    self.E.genotype_increase(gene_A, 1)
            else:
                self.B.genotype_increase(new_gene, 1)

#Production: Defines the first order production of one component by another component for example plasma cells producing antibodies.
class Production(OrderOne):
    def __init__(self, name, k, A, B):
        OrderOne.__init__(self, name, k, A)
        self.A = A
        self.B = B

    def react(self):
        base_gene = self.A.select_random()
        self.B.genotype_increase(base_gene, 1)

#Clearance: Defines the second order clearance of one component by another component through Ab/BCR binding. For example antibody-based clearance of an antigen
class Clearance(OrderTwoPhenotype):
    def __init__(self, name, k, A, B, max_rate, aff_factor, type):
        OrderTwoPhenotype.__init__(self, name, k, A, B, max_rate, aff_factor, type)
        self.A = A
        self.B = B

    def react(self):
        self.A.decrease(1)

#Replication: Defines the first order replication of one component into two identical components.
class Replication(OrderOne):
    def __init__(self, name, k, A):
        self.A = A
        OrderOne.__init__(self, name, k, A)

    def react(self):
        self.A.increase(1)

####################
#SYSTEM FUNCTIONS. #
####################
#TotalReaction: Contains the entire set of immune reactions that define the system. Calculates the reaction rate and Monte Carlo time step based on the Gillespie algorithm
class TotalReaction:
    def __init__(self):
        self.reaction_list = []
        self.n = 0
        self.__rate = []
        self.__total_rate = float(0.0)

    def add_reaction(self, reaction):
        self.reaction_list.append(reaction)
        self.n = len(self.reaction_list)

    def rate(self):
        self.__rate = []
        for i in range(0, self.n):
            self.__rate.append(self.reaction_list[i].rate())
        self.__total_rate = sum(self.__rate)

    def MC_TimeStep(self):
        r = random.random()
        self.rate()
        dt = (1 / self.__total_rate) * math.log(1 / r)
        #dt = float(2.0/1000)
        return dt

    def MC_React(self):
        r = random.random()
        outcome = float(0.0)
        for i in range(0, self.n):
            outcome += (self.__rate[i] / self.__total_rate)
            if (r <= outcome):
                self.reaction_list[i].react()
                break

#Vaccine
class Vaccine:
	def __init__(self):
		self.PopulationList = []
	
	def increase_ag1(self, ag1_inital, ag_range):
		increase_value = ag1_inital / 12
		for i in range(ag_range):
			self.data.append(increase_value)

#FileOutput: Outputs the simulation data into a data file, at defined time intervals. The simulation data is mainly the individual populations of each component of the system, 
# as well as subpopulations broken down by phenotype. 
class FileOutput:
    def __init__(self, data_file, step, population, antigen_list):
        """ :type self: object """
        self.time = float(0.0)
        self.data_file = data_file
        self.step = step
        self.population = population
        self._antigen_list = antigen_list
        self._antigen1 = self._antigen_list[0]#Cal09
        self._antigen2 = self._antigen_list[1]#BR07
        self._antigen3 = self._antigen_list[2]#SC18
        self._antigen4 = self._antigen_list[3]#SI06
        self._antigen5 = self._antigen_list[4]#test
        self._antigen6 = self._antigen_list[5]#test
        self._antigen7 = self._antigen_list[6]#test
        self._antigen8 = self._antigen_list[7]#test
        self._antigen9 = self._antigen_list[8]#test
        self._antigen10 = self._antigen_list[9]#test
        self._antigen11 = self._antigen_list[10]#test
        self._antigen12 = self._antigen_list[11]#test

    def write(self, time):
        if (time == 0.0 or time > (self.time + self.step)):
            line_out = ''
            line_out += str(time) + ';'
            line_out += 'V;'
            for i in range(7, len(self.population)):
                line_out += str(self.population[i].size()) + ';'
            line_out += 'B;'
            line_out += str(self.population[0].size()) + ';'  #GC
            line_out += str(self.population[1].size()) + ';'  #Stim
            line_out += str(self.population[2].size()) + ';'  #Memory
            line_out += str(self.population[3].size()) + ';'  #Plasma
            line_out += str(self.population[4].size()) + ';'  #Antibody
            line_out += str(self.population[5].size()) + ';'  #Naive
            line_out += str(self.population[6].size()) + ';'  #LLpC

            line_out += 'Bphen;'
            for i in range(1, 8):
                output = []
                output.append(self.population[0].phenotype_size(i, self._antigen1.ID()))
                output.append(self.population[1].phenotype_size(i, self._antigen1.ID()))
                output.append(self.population[2].phenotype_size(i, self._antigen1.ID()))
                line_out += str(max(output)) + ';'

            line_out += 'ABphen;'
            for i in range(1, 8):
                line_out += str(self.population[4].phenotype_size(i, self._antigen1.ID())) + ';'

            line_out += 'Nepi;'
            for i in range(0, 10):#updated 4/18/2023 to output epitopes in AG2
                line_out += str(self.population[5].epitope_size(i)) + ';'
            line_out += 'Bepi;'
            for i in range(0, 10):
                line_out += str(self.population[1].epitope_size(i)) + ';'
            line_out += 'Mepi;'
            for i in range(0, 10):
                line_out += str(self.population[2].epitope_size(i)) + ';'
            line_out += 'ABepi;'
            for i in range(0, 10):
                line_out += str(self.population[4].epitope_size(i)) + ';'

            line_out += 'ABcross_A1_A1;'
            output = self.population[4].calc_crossreactivity(self._antigen1, self._antigen1)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'
            line_out += 'ABcross_A2_A2;'
            output = self.population[4].calc_crossreactivity(self._antigen2, self._antigen2)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'
            line_out += 'ABcross_A3_A3;'
            output = self.population[4].calc_crossreactivity(self._antigen3, self._antigen3)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'
            line_out += 'ABcross_A4_A4;'
            output = self.population[4].calc_crossreactivity(self._antigen4, self._antigen4)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'
            line_out += 'ABcross_A5_A5;'
            output = self.population[4].calc_crossreactivity(self._antigen5, self._antigen5)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'
            line_out += 'ABcross_A6_A6;'
            output = self.population[4].calc_crossreactivity(self._antigen6, self._antigen6)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'
            line_out += 'ABcross_A7_A7;'
            output = self.population[4].calc_crossreactivity(self._antigen7, self._antigen7)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'
            line_out += 'ABcross_A8_A8;'
            output = self.population[4].calc_crossreactivity(self._antigen8, self._antigen8)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'
            line_out += 'ABcross_A9_A9;'
            output = self.population[4].calc_crossreactivity(self._antigen9, self._antigen9)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'
            line_out += 'ABcross_A10_A10;'
            output = self.population[4].calc_crossreactivity(self._antigen10, self._antigen10)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'
            line_out += 'ABcross_A11_A11;'
            output = self.population[4].calc_crossreactivity(self._antigen11, self._antigen11)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'
            line_out += 'ABcross_A12_A12;'
            output = self.population[4].calc_crossreactivity(self._antigen12, self._antigen12)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'
            line_out += 'ABcross_A1_A2;'
            output = self.population[4].calc_crossreactivity(self._antigen1, self._antigen2)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'
            line_out += 'ABcross_A1_A3;'
            output = self.population[4].calc_crossreactivity(self._antigen1, self._antigen3)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'
            line_out += 'ABcross_A1_A4;'
            output = self.population[4].calc_crossreactivity(self._antigen1, self._antigen4)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'
            line_out += 'ABcross_A1_A5;'
            output = self.population[4].calc_crossreactivity(self._antigen1, self._antigen5)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'
            line_out += 'ABcross_A1_A6;'
            output = self.population[4].calc_crossreactivity(self._antigen1, self._antigen6)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'
            line_out += 'ABcross_A1_A7;'
            output = self.population[4].calc_crossreactivity(self._antigen1, self._antigen7)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'
            line_out += 'ABcross_A1_A8;'
            output = self.population[4].calc_crossreactivity(self._antigen1, self._antigen8)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'
            line_out += 'ABcross_A1_A9;'
            output = self.population[4].calc_crossreactivity(self._antigen1, self._antigen9)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'
            line_out += 'ABcross_A1_A10;'
            output = self.population[4].calc_crossreactivity(self._antigen1, self._antigen10)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'
            line_out += 'ABcross_A1_A11;'
            output = self.population[4].calc_crossreactivity(self._antigen1, self._antigen11)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'
            line_out += 'ABcross_A1_A12;'
            output = self.population[4].calc_crossreactivity(self._antigen1, self._antigen12)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'
            line_out += str(output[2]) + ';'

            line_out += 'ABtrans;'
            output = self.population[4].calc_transcend(self._antigen_list)
            #for i in range(0, len(self._antigen_list)+1):
            for i in range(0, len(self._antigen_list)):
                line_out += str(output[i]) + ';'

            line_out += 'Btrans;'
            output0 = self.population[0].calc_transcend(self._antigen_list)
            output1 = self.population[1].calc_transcend(self._antigen_list)
            output2 = self.population[2].calc_transcend(self._antigen_list)
            #for i in range(0, len(self._antigen_list)+1):
            for i in range(0, len(self._antigen_list)):
                output = []
                output.append(output0[i])
                output.append(output1[i])
                output.append(output2[i])
                line_out += str(max(output)) + ';'

            line_out += 'ABneut;'
            output = self.population[4].calc_neutralization(self._antigen1, self._antigen2)
            line_out += str(output[0]) + ';'
            line_out += str(output[1]) + ';'

            line_out += 'ABdep;'
            output = self.population[4].calc_depletion(self._antigen1, self._antigen2)
            line_out += str(output) + ';'

            line_out += 'Bdiv;'
            output = []
            n = 0.25
            output.append(self.population[0].diversity2(n))
            output.append(self.population[1].diversity2(n))
            output.append(self.population[2].diversity2(n))
            line_out += str(max(output)) + ';'

            n = 0.50
            output.append(self.population[0].diversity2(n))
            output.append(self.population[1].diversity2(n))
            output.append(self.population[2].diversity2(n))
            line_out += str(max(output)) + ';'

            n = 0.75
            output.append(self.population[0].diversity2(n))
            output.append(self.population[1].diversity2(n))
            output.append(self.population[2].diversity2(n))
            line_out += str(max(output)) + ';'

            n = 1.0
            output.append(self.population[0].diversity2(n))
            output.append(self.population[1].diversity2(n))
            output.append(self.population[2].diversity2(n))
            line_out += str(max(output)) + ';'

            line_out += 'AB_ep1_phen_Vall;'
            output = self.population[4].calc_cross(self._antigen1, 0)
            for i in range(0, 4):
                line_out += str(output[i]) + ';'

            line_out += 'AB_ep2_phen_V0;'
            output = self.population[4].calc_cross(self._antigen1, 1)
            for i in range(0, 4):
                line_out += str(output[i]) + ';'

            line_out += 'AB_ep2_phen_V1;'
            output = self.population[4].calc_cross(self._antigen2, 1)
            for i in range(0, 4):
                line_out += str(output[i]) + ';'

            line_out += 'AB_ep2_phen_V2;'
            output = self.population[4].calc_cross(self._antigen3, 1)
            for i in range(0, 4):
                line_out += str(output[i]) + ';'

            line_out += 'AB_ep2_phen_V3;'
            output = self.population[4].calc_cross(self._antigen4, 1)
            for i in range(0, 4):
                line_out += str(output[i]) + ';'

            line_out += 'AB_ep3_phen_V0;'
            output = self.population[4].calc_cross(self._antigen1, 2)
            for i in range(0, 4):
                line_out += str(output[i]) + ';'

            line_out += 'AB_ep3_phen_V1;'
            output = self.population[4].calc_cross(self._antigen2, 2)
            for i in range(0, 4):
                line_out += str(output[i]) + ';'

            line_out += 'AB_ep3_phen_V2;'
            output = self.population[4].calc_cross(self._antigen3, 2)
            for i in range(0, 4):
                line_out += str(output[i]) + ';'

            line_out += 'AB_ep3_phen_V3;'
            output = self.population[4].calc_cross(self._antigen4, 2)
            for i in range(0, 4):
                line_out += str(output[i]) + ';'

            line_out += 'AB_ep4_phen_V0;'
            output = self.population[4].calc_cross(self._antigen1, 3)
            for i in range(0, 4):
                line_out += str(output[i]) + ';'

            line_out += 'AB_ep4_phen_V1;'
            output = self.population[4].calc_cross(self._antigen2, 3)
            for i in range(0, 4):
                line_out += str(output[i]) + ';'

            line_out += 'AB_ep4_phen_V2;'
            output = self.population[4].calc_cross(self._antigen3, 3)
            for i in range(0, 4):
                line_out += str(output[i]) + ';'

            line_out += 'AB_ep4_phen_V3;'
            output = self.population[4].calc_cross(self._antigen4, 3)
            for i in range(0, 4):
                line_out += str(output[i]) + ';'

            line_out += 'AB_ep5_phen_V0;'
            output = self.population[4].calc_cross(self._antigen1, 4)
            for i in range(0, 4):
                line_out += str(output[i]) + ';'

            line_out += 'AB_ep5_phen_V1;'
            output = self.population[4].calc_cross(self._antigen2, 4)
            for i in range(0, 4):
                line_out += str(output[i]) + ';'

            line_out += 'AB_ep5_phen_V2;'
            output = self.population[4].calc_cross(self._antigen3, 4)
            for i in range(0, 4):
                line_out += str(output[i]) + ';'

            line_out += 'AB_ep5_phen_V3;'
            output = self.population[4].calc_cross(self._antigen4, 4)
            for i in range(0, 4):
                line_out += str(output[i]) + ';'

            line_out += 'AB_ep6_phen_V0;'
            output = self.population[4].calc_cross(self._antigen1, 5)
            for i in range(0, 4):
                line_out += str(output[i]) + ';'

            line_out += 'AB_ep6_phen_V1;'
            output = self.population[4].calc_cross(self._antigen2, 5)
            for i in range(0, 4):
                line_out += str(output[i]) + ';'

            line_out += 'AB_ep6_phen_V2;'
            output = self.population[4].calc_cross(self._antigen3, 5)
            for i in range(0, 4):
                line_out += str(output[i]) + ';'

            line_out += 'AB_ep6_phen_V3;'
            output = self.population[4].calc_cross(self._antigen4, 5)
            for i in range(0, 4):
                line_out += str(output[i]) + ';'

            print(line_out)
            line_out += '\n'
            self.f.write(line_out)
            self.time = time

    def start(self):
        self.f = open(self.data_file, 'w')
        line_out = ''
        line_out += "Time;"

        line_out += 'V;'
        for i in range(7, len(self.population)):
            line_out += str(self.population[i].name()) + ';'

        line_out += 'B;'
        line_out += str(self.population[0].name()) + ';'
        line_out += str(self.population[1].name()) + ';'
        line_out += str(self.population[2].name()) + ';'
        line_out += str(self.population[3].name()) + ';'
        line_out += str(self.population[4].name()) + ';'
        line_out += str(self.population[5].name()) + ';'
        line_out += str(self.population[6].name()) + ';'
        line_out += 'Bphen;'

        tag = "B_stim"
        for i in range(1, 8):
            line_out += tag + "_" + str(i) + ';'
        #tag = "B_mem"
        #for i in range(1,8):
        #	line_out += tag+"_"+str(i)+';'
        line_out += 'ABphen;'
        tag = "AB"
        for i in range(1, 8):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'Nepi;'
        tag = "B_naive_AG"
        for i in range(0, 10): #updated 4/23/2023 to deal with AG2
            line_out += tag + "_" + str(i) + ';'
        line_out += 'Bepi;'
        tag = "B_stim_AG"
        for i in range(0, 10):#updated 4/23/2023 to deal with AG2
            line_out += tag + "_" + str(i) + ';'
        line_out += 'Mepi;'
        tag = "B_mem_AG"
        for i in range(0, 10):#updated 4/23/2023 to deal with AG2
            line_out += tag + "_" + str(i) + ';'
        line_out += 'ABepi;'
        tag = "AB_AG"
        for i in range(0, 10):#updated 4/23/2023 to deal with AG2
            line_out += tag + "_" + str(i) + ';'

        line_out += 'ABcross_A1_A1;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'
        line_out += 'ABcross_A2_A2;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'
        line_out += 'ABcross_A3_A3;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'
        line_out += 'ABcross_A4_A4;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'
        line_out += 'ABcross_A5_A5;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'
        line_out += 'ABcross_A6_A6;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'
        line_out += 'ABcross_A7_A7;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'
        line_out += 'ABcross_A8_A8;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'
        line_out += 'ABcross_A9_A9;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'
        line_out += 'ABcross_A10_A10;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'
        line_out += 'ABcross_A11_A11;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'
        line_out += 'ABcross_A12_A12;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'
        line_out += 'ABcross_A1_A2;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'
        line_out += 'ABcross_A1_A3;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'
        line_out += 'ABcross_A1_A4;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'
        line_out += 'ABcross_A1_A5;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'
        line_out += 'ABcross_A1_A6;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'
        line_out += 'ABcross_A1_A7;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'
        line_out += 'ABcross_A1_A8;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'
        line_out += 'ABcross_A1_A9;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'
        line_out += 'ABcross_A1_A10;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'
        line_out += 'ABcross_A1_A11;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'
        line_out += 'ABcross_A1_A12;'
        tag = "AB_cross"
        line_out += tag + ';'
        tag = "AB_spec1"
        line_out += tag + ';'
        tag = "AB_spec2"
        line_out += tag + ';'

        line_out += 'ABtrans;'
        tag = "AB_trans"
        #for i in range(0, len(self._antigen_list)+1):
        for i in range(0, len(self._antigen_list)):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'Btrans;'
        tag = "B_trans"
        #for i in range(0, len(self._antigen_list)+1):
        for i in range(0, len(self._antigen_list)):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'ABneut;'
        tag = "AB_neut1"
        line_out += tag + ';'
        tag = "AB_neut2"
        line_out += tag + ';'

        line_out += 'ABdep;'
        tag = "AB_dep1"
        line_out += tag + ';'

        line_out += 'Bdiv;'
        tag = "Bdiv_25"
        line_out += tag + ';'
        tag = "Bdiv_50"
        line_out += tag + ';'
        tag = "Bdiv_75"
        line_out += tag + ';'
        tag = "Bdiv_100"
        line_out += tag + ';'

        line_out += 'AB_ep1_phen_Vall;'
        tag = "ABE1All"
        for i in range(0, 4):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'AB_ep2_phen_V0;'
        tag = "ABE2V0"
        for i in range(0, 4):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'AB_ep2_phen_V1;'
        tag = "ABE2V1"
        for i in range(0, 4):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'AB_ep2_phen_V2;'
        tag = "ABE2V2"
        for i in range(0, 4):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'AB_ep2_phen_V3;'
        tag = "ABE2V3"
        for i in range(0, 4):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'AB_ep3_phen_V0;'
        tag = "ABE3V0"
        for i in range(0, 4):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'AB_ep3_phen_V1;'
        tag = "ABE3V1"
        for i in range(0, 4):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'AB_ep3_phen_V2;'
        tag = "ABE3V2"
        for i in range(0, 4):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'AB_ep3_phen_V3;'
        tag = "ABE3V3"
        for i in range(0, 4):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'AB_ep4_phen_V0;'
        tag = "ABE4V0"
        for i in range(0, 4):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'AB_ep4_phen_V1;'
        tag = "ABE4V1"
        for i in range(0, 4):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'AB_ep4_phen_V2;'
        tag = "ABE4V2"
        for i in range(0, 4):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'AB_ep4_phen_V3;'
        tag = "ABE4V3"
        for i in range(0, 4):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'AB_ep5_phen_V0;'
        tag = "ABE5V0"
        for i in range(0, 4):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'AB_ep5_phen_V1;'
        tag = "ABE5V1"
        for i in range(0, 4):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'AB_ep5_phen_V2;'
        tag = "ABE5V2"
        for i in range(0, 4):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'AB_ep5_phen_V3;'
        tag = "ABE5V3"
        for i in range(0, 4):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'AB_ep6_phen_V0;'
        tag = "ABE6V0"
        for i in range(0, 4):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'AB_ep6_phen_V1;'
        tag = "ABE6V1"
        for i in range(0, 4):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'AB_ep6_phen_V2;'
        tag = "ABE6V2"
        for i in range(0, 4):
            line_out += tag + "_" + str(i) + ';'

        line_out += 'AB_ep6_phen_V3;'
        tag = "ABE6V3"
        for i in range(0, 4):
            line_out += tag + "_" + str(i) + ';'

        line_out += '\n'
        self.f.write(line_out)  #outputs header
        self.write(0.0)  #writes starting conditions

    def finish(self):
        self.f.close()

    def set_time(t):
        self.time = t

    def gene_out(self):
        gene_file = self.data_file + ".gen"
        self.f = open(gene_file, 'w')

        gene_list = self.population[2].return_gene_list()
        size_list = self.population[2].return_size_list()
        for i in range(0, len(gene_list)):
            line_out = ''
            line_out += str(gene_list[i]) + ';'
            line_out += str(GeneEpitope(gene_list[i], self._antigen_list)) + ';'
            line_out += str(size_list[i]) + ';'
            line_out += '\n'
            print(line_out)
            self.f.write(line_out)  #outputs header
        self.f.close()

