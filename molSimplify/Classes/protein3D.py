#  @file protein3D.py
#  Defines protein3D class and contains useful manipulation/retrieval routines.
#
#  Written by HJK Group
#
#  Dpt of Chemical Engineering, MIT

# imports
from math import sqrt
import os, io
from molSimplify.Classes.AA3D import AA3D
from molSimplify.Classes.atom3D import atom3D
from molSimplify.Classes.helpers import pStripSpaces
import gzip
from itertools import chain
import urllib.request as urllib
import requests
from bs4 import BeautifulSoup
import pandas as pd
import string

# no GUI support for now

waters = {'AHOH', 'CHOH', 'HOH', 'EHOH', 'DHOH', 'BHOH'}

heme_type = {'AHEC', 'HEC', 'BHEC', 'ACOH', 'BCOH', 'SIR', 'CMO', 'HEM', 'B12',
             'BW9', 'CNC', 'BJ1S', 'PNI', 'COB', 'AJ1R', 'CBY', "CLA", "AHEM",
             "BHEM", "BCL", "ABCL", "BBCL", 'SRM', "F43", "HDD", 'HEB', 'BCB',
             'ZNH', 'AZNH', 'BZNH', 'HNI', 'FEC', 'AFEC', 'VOV', 'BVOV', 'DHE',
             '7HE', '6HE', 'MH0', 'DEU', "NTE", "HKL", 'BVQ', 'OBL', 'HEA',
             'AHEA', 'BHEA', 'HEV', 'ACMO', 'MI9', 'CL7', 'HAS', 'RUR', 'BCMO',
             'CHL', '83L', 'CHEM', 'DHEM', 'HCO', '89R', 'A89R', 'B89R', 'MNH',
             'HNN', 'WXP', 'VEA', 'COH', '4HE', '9QQ'}

others = {'ANCO', 'NCO', 'BNCO', 'NO', 'BNO', 'ANO', 'BEF', 'CEOH', '1PE',
          'P6G', 'CYN', 'TA6', 'MLI', '2PE', '15P', 'SO4', 'ACT', 'DMS', 'ZL4',
          'BALF', 'ALF', 'PRP', 'APRP', 'CAC', 'BGQ', 'DBH', 'PO4', 'POP',
          'GOL', 'PGE', 'PPV', 'BPPV', 'PG6', 'GXT', 'ZN8', 'APEG', 'BPEG',
          'PER', 'DVT', 'PG4', 'MDN', 'CGU', 'BGC', 'NO3', 'OXY', 'SIN', 'AKG',
          "F3S", "FES", "IMD", "8K2", "BE7", "ABE7", "TRS", "B22", "B17", 'PC',
          "RKL", "CUK", "CUA", "EQQ", "BTB", "GB1", "9E8", '9E5', "9DH", 'BU1',
          'EDO', 'GN1', 'QG1', "AQG1", "QG4", 'QFY', 'Q77', 'PX7', 'Q7A', 'OH',
          'BUB', 'PYR', 'Q71', "CL", 'EHM', 'LFC', 'Q6D', 'Q64', 'VO4', 'AOS',
          'FEL', 'RNS', 'KGY', 'FMT', 'MGF', 'TRS', "ZSP", "OPE", "IHP", '5MY',
          '724', 'ACY', 'Z59', "SCN", 'CHD', "FLC", 'FX7', 'PX7', "FT8", "1IU",
          'DA4', 'SF4', "WCC", '9NB', "FMT", 'APPV', "JNN", "AJNN", "BJNN",
          '8XZ', "PLM", "AV50", "BV50", 'V50', 'KIF', "ATRS", "BTRS", "PPY",
          "XC3", "P2H", "HOA", "ASY", 'GPP', "PEG", '1OU', "1OV", '19N', '1OS',
          '1OT', 'P5Y', "ABH", "8P8", "CO3", 'CLF', "HCA", 'NFU', 'NMY', "COA",
          "CIT", 'HAE', 'BB2', "FUC", '8L2', '8L5', "8LB", "8KN", "8KK", "8JQ",
          '8GJ', '8GK', "OXL", "HCQ", "HCH", 'JC2', '8SH', 'FCN', "1KM", 'TZD',
          '7MH', '7MF', '7ML', '7MK', "CGOL", "DMJ", "ADMJ", "BDMJ", "HS6",
          "PEO", 'S7B', 'HS7', 'SO3', "MBY", "AMBY", "DPO", "BDPO", "BHO",
          "ABHO", "BBHO", 'AGV', "CAC", "ACAC", "AGOL", "BGOL", "TPP", "B3N",
          'AB3N', 'M3S', "MPO", "AMPO", 'BIX', "ZDC", "FEO", "SHH", "FVZ",
          "TDP", 'K0I', "GAL", '2KT', "WBA", "BTB", "ABTB", "BBTB", "MVL",
          '6XO', 'FCO', "PO3", 'XXG', "QMS", "MYR", "THV", "THY", 'L88', "TDK",
          "TFD", "ATFD", "BTFD", "BFMT", 'O2Y', 'GVM', 'LTG', 'SOL', 'FUD',
          'TFG', 'AQA', "FSX", 'BFSX', 'ASF4', 'BSF4', 'TBV', 'BG6', 'MBD',
          'XCC', 'CUV', 'BCUV', 'AXCC', '2M8', 'VN4', '3PG', 'AF3', 'MAN',
          'AMAN', 'BMA', 'BBMA', 'LMS', 'SHA', 'ASO4', 'BSO4', 'PBC', 'APBC',
          'KMH', 'AAE', 'NAG', 'YIV', 'PVP', 'BACT', '7GR', 'XI7', 'ETX',
          '79E', '7K0', 'A7MK', 'HJ5', '24F', '6QC', '2FP', 'HS5', 'HS1', 'ZN',
          'HS3', '6H0', '8JH', 'FZ1', '79F', 'THW', 'IPT', 'SLB', 'SPV', 'PGA',
          '6NN', '10B', '7TR', '0XW', 'KG7', '0ZD', 'PEP', 'APEP', 'BPEP',
          '2PG', 'B2PG', '6ML', 'B6ML', 'H2S', 'CH2S', 'AH2S', '0XX', '7AF',
          'SPV', 'AZI', 'NMN', 'W45', '4NC', 'DHY', 'NG1', 'IWX', '6EZ', 'FBP',
          'ACA', 'AACA', 'BACA', 'Y08', 'Y16', 'YZ6', 'Y10', '6C1', '2PN',
          'ARS', 'DPJ', 'OLC', 'PE5', 'COM', '294', 'IMH', 'AIMH', 'BIMH',
          'FOM', 'CO2', 'BZS', '7FY', 'B7FY', 'SRO', 'HTD', 'HTL', 'XYP',
          'ZLM', 'ZLS', 'ZLV', 'V1P', 'PGH', 'P33', 'OLW', 'K9L', 'L5Z', '6J0',
          'C65', 'HA6', '9RB', 'N5J', '7H1', 'CDPO', 'N4R', 'TSN', 'E1Z', 'MG',
          '6DK', 'YED', '5OO', "TLA", "QPT", "CPT", "HGB", "IL5", "WO6", "MOS",
          "MTE", '0F3', 'JSC', "BJSC", "JSD", "JSE", "AJSE", '3WB', "RKM",
          'HB8', '5YA', 'MMZ', 'XCZ', '82N', 'FE2', 'OXD', 'OE2', '5SR', 'GRX',
          'AIV', '0UE', 'MCO', 'X8Z', '5V3', '5V0', 'T34', 'GRZ', 'T86', 'GT2',
          'GRK', 'GQZ', 'GM5', '5PU', 'M2A', 'TG4', '5RD', 'TG5', 'ATG5',
          'BTG5', 'TD6', '5ON', '9FK', '5O5', '5OM', '5OG', '3YP', 'IP8', 'NI',
          'TPS', 'TMG', '55L', 'A55L', 'B55L', 'LTQ', '5LD', 'TD4', 'ATD4',
          'BTD4', 'E7Z', 'PPK', 'APPK', 'E7Q', '5Q1', 'O84', 'NO2', 'BRJ', 'F',
          'D4B', 'ZAR', 'AZAR', 'BZAR', 'XCH', 'D0Z', 'CXH', 'FZZ', 'AFZZ', 
          'BYN', 'BBYN', '4LU', 'A4LU', '4MJ', 'BX5', 'B9Z', 'OLC', 'LPW', 'O',
          'BJ5', 'BJ2', 'B9Z', 'VOH', 'B9N', 'J6A', 'J4P', 'J4S', 'J4V', 'J1V',
          'J2D', 'J0Y', 'J1G', 'WNN', 'HBV', 'HBJ', 'AHBJ', 'BHBJ', '7MT',
          'MMA', 'AMMA', 'BMMA', 'BNO2', 'MLT', 'BME', '8YL', 'A8YL', '8F2',
          '8F3', 'CAP', 'HZK', 'HZW', 'HZE', 'XUL', 'WS2', 'AWS2', '4BW',
          'EOJ', 'TQJ', 'FBM', 'FBJ', 'TP9', '1GP', 'BU3', 'AEDO', 'CQB', 'HG',
          'AAKG', '8GH', '8HK', '8HE', '8GE', 'BTN', 'V13', 'AV13', 'BV13',
          'S3C', 'BCN', 'MPD', 'A2G', 'PMV', 'ETA', 'BETA', 'FB2', '8Q8',
          'JTY', 'NHC', 'ANHC', 'BNHC', 'ENL', 'AENL', 'BENL', 'BES', '2NO',
          'VNI', '476', 'MYA', 'HZT', 'HZE', 'PEB', 'URE', '93W', '8V5', '3PO',
          'J2H', 'SOR', 'XLS', 'XYL', 'V14', 'EZL', 'DJ3', 'ADJ3', 'BDJ3',
          '91R', 'B91R', '91T', 'A91T', 'BO3', '1R5', 'ZP1', '92K', 'QSP',
          '1VQ', '9GB', '9GE', 'DHB', '1SD', 'TOR', 'E1E', 'E1F', 'AE1F',
          'BE1F', '3S0', 'Z3K','NHW', 'E1G', 'EA3', 'JTW', 'AJTW', 'BJTW',
          'BPGE', '9KB', '9KE', 'BTW', 'GAE', 'LA6', 'BLA6', '0FV', 'A0FV',
          'ALA6', 'FPS', 'IPM', '13P', 'BCT', 'VNT', 'J9Y', 'AJ9Y', 'ZZ7',
          'ZTP', 'ASIN', 'BSIN', 'TAR', 'S9N', 'BS9N', 'PK8', 'BPO4', 'CFO',
          'NTM', 'PS7', 'AHD', 'G6P', 'EF1', 'AEF1', 'BEF1', '4NG', 'LG6',
          'LXP', 'ADA', 'AADA', 'BADA', 'CADA', 'GTR', 'BGTR', 'YOM', '4KN',
          'PBD', 'PH2', '4J9', 'PUD', 'TGD', 'ATGD', 'BGTD', 'HY0', "4KD",
          "4LC", '4KB', '4KC', "R3N", "QZH", "NHE", "RAM", "GST", "FUL", "FUC",
          "AFUL", "BFUC", "BDF", "FMN", "CEL", "NMH", 'RF7', 'RC5', '43G',
          '2HA', 'A2HA', 'G16', "GCD", "GMH", '293', 'PIM', "XSP", '0V5', 'CD',
          'APV', '23N', "GSH", "HSM", "HCS", "AZM", 'C2E', '3X1', '3ZP', 'MVI',
          'VD9', 'MIV', 'ZM3', '3TV', 'HXD', 'GLC', 'THJ', 'NOE', 'BTHJ', 'PI',
          '402', '40E', 'HIO', 'NDH', 'T5X', 'T6F', 'FDP', 'MVI', '03W', '3UG',
          'A03W', '3UF', 'A3UF', 'B3UF', 'BPRP', '3PY', 'A3PY', 'B3PY', 'KAN',
          'SIA', 'AR6', 'ICT', '9L3', 'DO1', '865', '3S6', 'HBU', 'V1S', '786',
          'GOA', 'R5X', 'APO4', 'L6G', '8JS', 'BSI', 'V1D', 'GOX', 'OEV', 'CU',
          'GIM', 'TGJ', 'FC4', 'AFC4', '6KE', 'B6KE', '6ED', '6E7', '6EP',
          '6EN', 'MQN', 'TZM', 'MMK', 'TQU', 'N9W', 'MB7', '9NK', '9C2', '8XW',
          '8V8', '7O8', '7NU', 'P6T', '1KH', 'A1KH', 'XPE', 'OJ4', 'ACO3',
          'BCO3', 'D7K', '4OP', '4PS', 'APER', 'BPER', 'AAZI', 'ANO2', '8M0',
          'E43', 'TBR', 'WO4', 'AACT', 'ZRW', 'REO', 'CU1', 'MO7', 'TEW',
          'TFM', '1CL', 'GXW', 'ZW2', '2HE', 'ARW', 'CFM', 'IKX', 'DMU', 'PEK',
          'TGL', 'BAZI', 'GDS', '9RU', 'A9RU', 'ACOM', '38I', 'ETS', 'A8XQ',
          '8XQ', 'DZZ', 'TZZ', '1GO', 'T4F', 'T5F', '7BH', '7BF', 'GVR', 'CYC',
          'BIE', 'AGLC', 'BBGC', 'ETF', 'AETF', 'BETF', '752', 'BLMR', 'ALMR',
          'LMR', 'B3P', 'PBE', 'NB4', 'NH3', 'ACO', 'DAO', '1EZ', 'PHB', '7DV',
          'BEZ', 'RC4', 'DTT', '6YQ', 'V51', '6YZ', 'V26', 'V49', '6Z9', '5EF',
          'A6Z9', 'B6Z9', 'A6N', '6YH', 'A6YH', 'B6YH', 'ICS', '2YU', '4FH',
          '5MH', '3MH', '3FH', 'ANV', '6XD', 'KKK', 'M8V', 'M8J', '642', 'X7A',
          'MFU', 'MO0', 'FRU', 'CV50', 'L8J', 'MZM', 'PYZ', 'DOR', '5LC', 'SE',
          'YNC', 'FUS', 'M3Q', 'KED', 'T6Z', 'YKG', 'XJE', 'FQV', 'FQY', 'VKE',
          'TWB', 'HL6', 'EZ1', '3ES', 'EPE', 'AHN', 'UCZ', '9P3', '9O9', '9OU',
          '9OC', 'M9F', 'E08', 'PNN', 'ZA3', 'FM5', '34F', 'PDO', 'E8S', 'ACE',
          'BUA', 'ABUA', 'BBUA', 'B13P', 'NY2', 'CL6', 'WPG', 'IMP', 'DQJ',
          'DQG', '90V', 'DNV', 'DKV', 'DKS', 'DKY', 'DKP', 'X93', 'BNAG', 'MO',
          'NDG', 'ANDG', 'BNF', 'FXY', 'TPU', 'NIO', 'VGL', 'WQP', 'A40E',
          'BHIO', 'II4', 'NCA', 'SV2', 'CQM', 'AIMD', 'BIMD', '3L4', 'FBV',
          'C6J', 'C6M', 'FBW', 'IUG', 'RMN', 'OKG', 'VVO', 'NH4', 'MOO', 'CWO',
          'PLL', 'AMOO', 'LJB', '47S', 'A47S', 'RRE', 'QHL', 'AQHL', 'R1C',
          'NFS', 'DMA', 'PAY', 'PQQ', 'SGN', 'IDS', 'GCU', 'BEDO', 'R7U', 'RU',
          'R9A', '11R', 'R5A', 'AR5A', 'R5B', 'BR5B', 'R6A', 'R4A', 'AR4A',
          'BR4A', 'M27', '6LL', 'PCD', 'H1Z', 'ICE', '4WA', 'TPQ', 'ATPQ',
          'A6R', 'LHW', 'M10', 'ANO3', 'PDC', 'APDC', 'BPDC', 'OMO', 'BEMC',
          'EMC', 'DO3', 'FC6', 'AJSD', 'AJSC', 'GXZ', 'AFBU', 'BFBU', 'FBU',
          'C2O', '28T', 'A28T', 'B28T', 'KJS', 'RZB', 'RZE', 'TAM', 'CUO',
          '4JE', 'S20', 'AS20', 'BS20', '4J0', '4JC', '1PT', '0FZ', 'D7A',
          'B09', 'A09', 'RU1', 'SMO', '9UX', 'UNX', 'RCS', 'ARCS', 'SRX',
          'MMC', 'APEO', 'BPEO', 'TPT', 'EMT', 'CU6', 'MSS', '6RP', 'EDT',
          'KLT', 'MSS', '9ZQ', 'R1Z', 'BQHL', 'ORS', 'SI7', 'SI8', 'ASQ1',
          'DPY', 'BSQ1', 'SQ1', 'DRP', '4IR', 'S32', 'S18', 'IUM', 'PMB', 'FE',
          'ACPT', 'RHD', 'IRI', 'AST', '31Q', 'RUU', 'C1O', 'DCD', 'QUE', 'MN'
          'KMP', 'A1PT', 'IOD', 'BIOD', 'AQPT', 'ACPT', 'BCPT', '6BP', 'ACL',
          'B9RU', '3V9', 'HM6', 'GHE', 'F6Q', 'AF6Q', 'BF6Q', 'BAST', 'PSC',
          'PGV', 'CPO4', 'CDL', 'APCD', '2MO', 'EFK', 'MM2', 'A12', 'AJCT',
          '52G', 'A52G', 'B52G', 'JCT', 'BJCT', 'UNL', 'DVW', 'ATPP', 'B8EL',
          '8EL', 'P23', '8FL', '8PO', '8EF', '8N9', '4LF', 'PT7', 'BNZ', 'ART',
          'AART', 'BART', 'C2C', '7G4', '7CZ', '0TR', 'M3T', 'P9G', 'AP9G',
          'BP9G', '4H2', '4HO', 'HQT', '7FH', 'B7FH', 'HT4', '75Y', '75W',
          '2UJ', '108', '109', 'BEW', '1SA', 'BO1', 'BOS', '2T8', 'PBQ', 'KCS',
          'A2T8', 'B19', 'SBW', 'YPT', 'AYPT', 'BYPT', 'TM7', 'PO1', '2J0',
          '25Y', 'RU7', '498', 'J3K', 'L7T', 'L4K', 'L1Q', 'KJ5', 'KFH', 'KGK',
          'AKGK', 'BKGK', 'KBZ', 'FBT', 'AFBT', 'BFBT', 'KBW', 'KBB', 'SAN',
          '45L', 'FBS', 'K5W', '4NZ', '95B', '9RT', 'RIR', '9QB', 'RU2', '9Q8',
          'A9Q8', 'FL1', 'V21', 'JR3', 'ACU1', 'BCU1', 'CACY', 'CDJ3', 'AA6R',
          'J0N', 'JR8', 'J0K', 'REQ', 'R1N', 'O1N', 'OT1', '0TN', 'RKP', '58Z',
          'RCZ', 'ARCZ', 'BRCZ', 'LC1', 'B8EF', '8EO', 'CAD', 'WO5', 'A73M',
          '73M', 'CCPT', 'MTV', '3BS', 'COX', '6M4', 'A6M4', 'B6M4', '0WZ',
          '2OP', '6B7', 'AC1', 'BL0', 'KEL', 'RML', 'TUO', '0KA', 'A2C', 'GZC',
          'GUS', '5WM', '5WN', 'H4O', 'DOD', 'GP9', 'YX0', 'WR2', 'CM9', 'CMW',
          'CII', 'ACM9', 'BCM9', 'ACII', 'BCII', 'TDC', 'DQY', 'M11', 'Q39',
          '1BN', 'FZ0', 'FYL', 'HHR', 'BPY', 'ZI6', 'ZI7', 'XBP', 'FQR', 'FQU',
          'A13P', 'LGU', 'FPY', 'CAQ', 'M2C', 'M3C', 'BFQ', '1AE', 'M1C', 'CS',
          'FMF', 'AFMF', 'BFMF','GUL', 'DFS', 'DFB', '164', '24B', 'LKA', 'NA',
          'LKS', '4SN', '4TE', '5L2', 'F4S', 'NFV', 'S0N', 'OEX', 'LMG', 'KC6',
          'P6F', 'ATSN', 'BTSN', 'PPF', 'HIW', '210', 'E90', 'E65', 'E2I',
          'E50', 'BVA', 'RZ1', 'RZ0', 'NK1', 'NK2', 'HPY', 'QKK', 'LB1', 'AUD',
          '3D1', '86B', 'RYZ', 'RYY', 'RYV', 'RZ8', 'RZ7', 'KPL', 'NLA', 'E49',
          'E59', 'AE59', 'BE59', 'EVJ', 'INS', 'AINS', 'BINS', '144', 'TCW',
          'EVI', 'EVH', 'EVG', 'EVF', 'BM3', 'IPD', 'EVE', 'EVD', 'LAC', 'FUQ',
          '4J8', 'F3K', 'AF3K', 'BF3K', 'ICG', 'CSO4', 'ZKG', '8RO', 'CDO3',
          'ADO3', 'BDO3', 'WO3', 'CRG', 'RHX', 'MOM', 'J9H', 'AEMC', '6U4',
          'QID', '42G', 'QEQ', 'PCR', 'HDE', 'BHDE', 'ORT', 'CZM', 'AHDD',
          'HDD', 'BHDD', 'CHDD', 'DHDD', 'BAKG', 'KYT', 'PRT', 'UZY', 'ARA',
          'C90', 'UBG', 'GLV', 'O81', 'ZOL', 'RZP', 'I8P', 'IBO', 'AIBO', 'BA',
          'BIBO', 'AFES', 'AGO', 'PSJ', '3E9', 'RHOH', 'NHOH', 'S3F', 'XS6',
          'SF3', 'XE1', 'AF3S', 'AFCO', '6YP', '3W6', '3W8', 'M4S', 'DCJ',
          'GWW', 'GWN', 'MM1', '16G', 'AICS', 'ICH', 'BICH', 'BGW', 'C4R',
          'BUNX', 'E3D', 'ACU', 'RUI', '4KV', 'BCU', '0JC', 'A0JC', 'B0JC',
          'KYS', 'AKYS', '8CY', 'A8CY', 'B8CY', 'LSI', 'CSF3', 'BSF3', 'BNO3',
          'DNC', 'BFES', 'H5Y', 'DP6', 'NGA', 'YOK', '5RP', 'GOS', '9TY',
          'PW8', '9YD', 'Z82', '0LX', 'PIS', '60N', '0S5', "IPR", '0WO', '0RV',
          'EF9', '9HK', 'GB3', "ON1", '84A', 'OY5', "BOY5", 'BH2S', 'ENO',
          'H1W', "AH1W", "BH1W", 'H1T', 'DHJ', "ETK", "IKY", 'U2P', 'VGS',
          'ARZ', 'V90', 'AV90', 'BV90'}

nucleotides = {'UTP', 'DUP', 'ADUP', 'BDUP', 'TXE', 'UNP', 'GDP', 'AMP', 'A66',
               '7DD', 'UPG', 'TNM', 'NAI', 'ANAI', 'ACP', 'I2C', 'DGI', 'DUD',
               'IDP', 'BA3', 'DGT', 'DAT', 'C38', 'UDP', 'NAD', 'B4P', 'UFT',
               'AGS', 'DUT', 'AU1', 'ADP', 'APR', 'GTP', '01G', 'ANP', 'B6P',
               'GMV', 'TTP', 'CTP', 'NCC', 'AP5', '8PA', 'APC', 'GPN', 'G5P',
               'UMP', 'ATP', 'TNV', "MGD", "MBO", "NYN", "BNYN", 'C5P', 'FAD',
               '0KX', 'DZ4', 'XG4' , 'D3T', 'GNP', 'Y9Z', 'GCP', "DTP", 'JSQ',
               "ADCP", "GH3", "G", 'E64', "ADTP", "ND7", "AP2", 'ZAN', "STP",
               "UD1", 'FE9', "ASTP", "BSTP", '6X9', 'T61', 'GQW', 'XG4', 'UDP',
               "AXG4", "BXG4", 'AS', 'A', '8DG', 'DAU', 'EO7', 'H84', 'BCTP',
               '2TM', 'DCT', 'U5P', 'CDP', '2KH', 'H27', '3BI', 'BHX', 'G1C',
               'DUN', 'DCM', 'A3P', 'GDU', 'HF7', 'HDV', 'WS1', 'GAV', 'DMV',
               '4S9', '4SZ', '4SY', '4TL', '4TM', 'G3D', 'DLL', 'VT6', '1FZ',
               'ATTP', 'AUDP', 'BUDP', "DCP", 'CA0', "UMA", '48O', "GSP",
               '4GW', 'B4GW', '3GW', "A3GW", '5GW', '2GW', '3AT', 'THP', 'ZDA',
               'AZDA', 'BZDA', 'DG3', 'DDS', '5GP', 'WS3', 'U37', '5AD', 'ODP',
               'THU', 'AMGD', 'BMGD', 'AMBO', 'BMBO', 'MCN', 'PGD', 'MD1',
               '5FA', 'AHZ', 'BTTP', 'QST', 'XJS', 'TMP', 'ATMP', 'AGTP',
               'BGTP', 'TAT', 'G2P', "BDCP", "H1Q", 'AU5P', 'BU5P'}

nonstandard_aas = {'ACAS', 'BCAS', 'CAS', 'TYD', 'MSO', 'MSE', 'CAF', 'HIS',
                   'AHO', 'BH2', 'NLG', 'MHS', "OGA", "CMH", "ACMH", "BCMH",
                   "NGH", 'CSO', 'INN', 'KCX', "ORN", 'CSX', "DGL", "DAS",
                   'BFD', "GLY", '6MG', 'CSS', '8OG', "16E", "ALY", '5H8',
                   "TOX", "ATOX", 'S17', 'SEP', "OCS", 'PSW', 'APSW', 'BPSW',
                   'E41', 'ARK', 'OKD', 'AR8', 'CI9', 'LEU', '1U4', '6NG',
                   'R47', 'Y38', '0PJ', 'OHN', 'AR47', 'BR47', 'R45', 'R4B',
                   'R4C', 'CSD', 'ACSX', 'AOCS', 'CSD', 'BL2', 'LEH', 'E37',
                   'ILE', 'VAL', 'APHE', 'PHE', 'TYR', 'TRP', 'TYQ', 'BTYQ',
                   'BC8', '5OL', '5H9', '56O', 'A56O', 'B56O', 'D0W', 'TPO',
                   '999', '4SQ', 'FGY', 'ACSD', 'HL5', '0Z7', 'B0Z7', 'DNP',
                   'BOCS', 'SHT', 'BSHT', 'C6L', 'TX4', 'CYS', '0X9', 'ASV',
                   "ACYS", 'BCYS', 'ALA', 'NHK', 'BINN', 'ACV', 'AGLY', 'BGLY',
                   '1K4', 'AOGA', 'BOGA', 'HIC', 'AHIC', 'BHIC', '2CO', 'MD6',
                   'UN9', '1B3', 'X16', '409', 'HSL', 'NCD', 'BCDH', 'ACDH',
                   'CDH', 'EVV', 'EAL', 'BCAF', 'M2W', 'WT4', 'AWT4', 'BWT4',
                   'ZZU', 'DAH', 'BDAH', 'DTR', 'HYP', 'ATYQ', 'ASEP', 'BSEP',
                   'BIR', 'VAZ', 'VB1', 'V20', 'GLU', 'FB6', 'FB5', 'LIZ',
                   'CUS', 'CME', 'SAH', 'BCV', 'ACSO', 'BCSO', 'LYS', 'BLEU',
                   'MET', 'HS8', 'AHS8', '2LT', 'BMET', 'S23', 'SAM', "ASAM",
                   "BSAM", "ACSS", "BCSS", "SCSS", "BCSD"}

class protein3D:
    """Holds information about a protein, used to do manipulations.  Reads
    information from structure file (pdb, cif) or is directly built from
    molsimplify.
    
    """
    
    def __init__(self, pdbfile='undef'):
        # Number of amino acids
        self.naas = 0 
        # Number of atoms not part of proteins
        self.nhetatms = 0
        # Number of chains
        self.nchains = 0
        # Set of amino acids
        self.aas = set()
        # Dictionary of all atoms
        self.atoms = {}
        # Dictionary of heteroatoms
        self.hetatms = {}
        # Dictionary of chains
        self.chains = {}
        # Dictionary of missing atoms
        self.missing_atoms = {}
        # List of missing amino acids
        self.missing_aas = []
        # List of possible AA conformations 
        self.conf = []
        # R value
        self.R = -1
        # Rfree value
        self.Rfree = -1
        # PDB file (if applicable)
        self.pdbfile = pdbfile
        # Holder for metals
        self.metals = False
        # Bonds
        self.bonds = {}
        # Data completeness
        self.DataCompleteness = 0
        # RSRZ value
        self.RSRZ = 100
        # TwinL score
        self.TwinL = 0
        # TwinL^2 score
        self.TwinL2 = 0
    
    def setAAs(self, aas):
        """ Set amino acids of a protein3D class to different amino acids.
        Parameters
        ----------
            aas : list
                Contains AA3D amino acids
        """
        self.aas = aas
        self.naas = len(aas)

    def setAtoms(self, atoms):
        """ Set atom indices of a protein3D class to atoms.
        
        Parameters
        ----------
            atoms : dictionary
                Keyed by atom index
                Valued by atom3D atom that has that index
        """
        self.atoms = atoms

    def setHetatms(self, hetatms):
        """ Set heteroatoms of a protein3D class to different heteroatoms.
        
        Parameters
        ----------
            hetatms : dictionary
                Keyed by heteroatoms
                Valued by a list containing the overarching molecule and chain
        """
        self.hetatms = hetatms
        self.nhetatms = len(hetatms.values())

    def setChains(self, chains):
        """ Set chains of a protein3D class to different chains.

        Parameters
        ----------
            chains : dictionary
                Keyed by desired chain IDs.
                Valued by the list of AA3D amino acids in the chain.
        """
        self.chains = chains
        self.nchains = len(chains.keys())

    def setMissingAtoms(self, missing_atoms):
        """ Set missing atoms of a protein3D class to a new dictionary.

        Parameters
        ----------
            missing_atoms : dictionary
                Keyed by amino acid residues of origin
                Valued by missing atoms
        """
        self.missing_atoms = missing_atoms

    def setMissingAAs(self, missing_aas):
        """ Set missing amino acids of a protein3D class to a new list.

        Parameters
        ----------
            missing_aas : list
                List of missing amino acids.
        """
        self.missing_aas = missing_aas
            
    def setConf(self, conf):
        """ Set possible conformations of a protein3D class to a new list.

        Parameters
        ----------
            conf : list
                List of possible conformations for applicable amino acids.
        """
        self.conf = conf
            
    def setR(self, R):
        """ Set R value of protein3D class.

        Parameters
        ----------
            R : float
                The desired new R value.
        """
        self.R = R
            
    def setRfree(self, Rfree):
        """ Set Rfree value of protein3D class.

        Parameters
        ----------
            Rfree : float
                The desired new Rfree value.
        """
        self.Rfree = Rfree

    def setRSRZ(self, RSRZ):
        """ Set RSRZ score of protein3D class.

        Parameters
        ----------
            RSRZ : float
                The desired new RSRZ score.
        """
        self.RSRZ = RSRZ
            
    def getMissingAtoms(self):
        """ Get missing atoms of a protein3D class.

        """
        return self.missing_atoms.values()
    
    def getMissingAAs(self):
        """ Get missing amino acid residues of a protein3D class.

        """
        return self.missing_aas
    
    def countAAs(self):
        """ Return the number of amino acid residues in a protein3D class.

        """
        return self.naas

    def findAtom(self, sym="X"):
        """
        Find atoms with a specific symbol that are contained in amino acids.
        
        Parameters
        ----------
            sym: str
                element symbol, default as X.

        Returns
        ----------
            inds: list
                a list of atom index with the specified symbol.
        """
        inds = []
        for aa in self.aas:
            for (ii, atom) in aa.atoms:
                if atom.symbol() == sym:
                    inds.append(ii)
        return inds

    def findHetAtom(self, sym="X"):
        """
        Find heteroatoms with a specific symbol.
        
        Parameters
        ----------
            sym: str
                element symbol, default as X.

        Returns
        ----------
            inds: list
                a list of atom index with the specified symbol.
        """
        inds = []
        for ii, atom in enumerate(self.hetatms.keys()):
            if atom.symbol() == sym:
                inds.append(ii)
        return inds

    def findAA(self, three_lc="XAA"):
        """
        Find amino acids with a specific three-letter code.
        
        Parameters
        ----------
            three_lc: str
                three-letter code, default as XAA.

        Returns
        -------
            inds: set
                a set of amino acid indices with the specified symbol.
        """
        inds = set()
        for aa in self.aas:
            if aa.three_lc == three_lc:
                inds.add((aa.chain, aa.id))
        return inds

    def getChain(self, chain_id):
        """ Takes a chain of interest and turns it into its own protein3D class instance.

        Parameters
        ----------
            chain_id : string
                The letter name of the chain of interest

        Returns
        -------
            p : protein3D
                A protein3D instance consisting of just the chain of interest
        """
        p = protein3D()
        p.setChains({chain_id: self.chains[chain_id]})
        p.setAAs(set(self.chains[chain_id]))
        p.setR(self.R)
        p.setRfree(self.Rfree)
        missing_aas = []
        for aa in self.missing_aas:
            if aa.chain == chain_id:
                missing_aas.append(aa)
        p.setMissingAAs(missing_aas)
        gone_atoms = {}
        for aa in self.missing_atoms.keys():
            if aa.chain == chain_id:
                gone_atoms[aa] = self.missing_atoms[aa]
        p.setMissingAtoms(gone_atoms)
        gone_hets = self.hetatms
        atoms = {}
        for a_id in self.atoms:
            aa = self.getResidue(a_id)
            if aa != None:
                if aa.chain == chain_id:
                    atoms[a_id] = self.atoms[a_id]
            else:
                if chain_id not in gone_hets[(a_id, self.atoms[a_id])]:
                    del gone_hets[(a_id, self.atoms[a_id])]
                else:
                    atoms[a_id] = self.atoms[a_id]
        p.setHetatms(gone_hets)
        p.setAtoms(atoms)
        bonds = {}
        for a in self.bonds.keys():
            if a in p.atoms.values():
                bonds[a] = set()
                for b in self.bonds[a]:
                    if b in p.atoms.values():
                        bonds[a].add(b)
        p.setBonds(bonds)
        return p

    def getResidue(self, a_id):
        """ Finds the amino acid residue that the atom is contained in.

        Parameters
        ----------
            a_id : int
                the index of the desired atom whose residue we want to find

        Returns
        -------
            aa : AA3D
                the amino acid residue containing the atom
                returns None if there is no amino acid
        """
        for aa in self.aas:
            if (a_id, self.atoms[a_id]) in aa.atoms:
                return aa
        for aa in self.missing_atoms.keys():
            if (a_id, self.atoms[a_id]) in self.missing_atoms[aa]:
                return aa
        return None # the atom is a heteroatom

    def stripAtoms(self, atoms_stripped):
        """ Removes certain atoms from the protein3D class instance.
        
        Parameters
        ----------
            atoms_stripped : list
                list of atom3D indices that should be removed
        """
        for aa in self.aas:
            for (a_id, atom) in aa.atoms:
                if a_id in atoms_stripped:
                    self.aas[aa].remove((a_id, atom))
                    atoms_stripped.remove(a_id)
        for (h_id, hetatm) in self.hetatms.keys():
            if h_id in atoms_stripped:
                del self.hetatms[(h_id, hetatm)]   
                atoms_stripped.remove(h_id)

    def stripHetMol(self, hetmol):
        """ Removes all heteroatoms part of the specified heteromolecule from
            the protein3D class instance.

        Parameters
        ----------
            hetmol : str
                String representing the name of a heteromolecule whose
                heteroatoms should be stripped from the protein3D class instance
        """
        for (h_id, hetatm) in self.hetatms.keys():
            if hetmol in self.hetatms[(h_id, hetatm)]:
                del self.hetatms[(h_id, hetatm)] 

    def findMetal(self, transition_metals_only=True):
        """Find metal(s) in a protein3D class.
        Parameters
        ----------
            transition_metals_only : bool, optional
                Only find transition metals. Default is true.
                
        Returns
        -------
            metal_list : list
                List of indices of metal atoms in protein3D.
        """
        if not self.metals:
            metal_list = []
            for (i, atom) in self.hetatms.keys(): # no metals in AAs
                if atom.ismetal(transition_metals_only=transition_metals_only):
                    metal_list.append(i)
            self.metals = metal_list
        return (self.metals)

    def freezeatom(self, atomIdx):
        """Set the freeze attribute to be true for a given atom3D class.

        Parameters
        ----------
            atomIdx : int
                Index for atom to be frozen.
        """

        self.atoms[atomIdx].frozen = True

    def freezeatoms(self, Alist):
        """Set the freeze attribute to be true for a given set of atom3D classes, 
        given their indices. Preserves ordering, starts from largest index.

        Parameters
        ----------
            Alist : list
                List of indices for atom3D instances to remove.
        """

        for h in sorted(Alist, reverse=True):
            self.freezeatom(h)

    def getAtom(self, idx):
        """Get atom with a given index.

        Parameters
        ----------
            idx : int
                Index of desired atom.

        Returns
        -------
            atom : atom3D
                atom3D class for element at given index.

        """

        return self.atoms[idx]

    def getBoundAAs(self, h_id):
        """Get a list of amino acids bound to a heteroatom, usually a metal.

        Parameters
        ----------
            h_id : int
                the index of the desired (hetero)atom origin

        Returns
        -------
            bound_aas : list
                list of AA3D instances of amino acids bound to hetatm
        """
        bound_aas = []
        for b_id in self.atoms.keys():
            b = self.atoms[b_id]
            if self.atoms[h_id] not in self.bonds.keys():
                g = self.atoms[h_id].greek
                s = self.atoms[h_id].sym
                if g == 'A' + s or g == 'B' + s: # just different conformation
                    return None
            elif b in self.bonds[self.atoms[h_id]]:
                if self.getResidue(b_id) != None:
                    bound_aas.append(self.getResidue(b_id))
                elif (b_id, b) in self.hetatms.keys():
                    # accommodate nonstandard amino acids
                    b_mol = self.hetatms[(b_id, b)][0]
                    #if b_mol in nonstandard_aas:
                    #print(self.atoms[h_id].sym, b_mol)
                    if b_mol not in waters and b_mol not in nucleotides and b_mol not in heme_type and b_mol not in others:
                        if b_mol not in nonstandard_aas:
                            print('UNKNOWN LIGAND', b_mol)
        return bound_aas
    
    def readfrompdb(self, text):
        """ Read PDB into a protein3D class instance.

        Parameters
        ----------
            text : str
                String of path to PDB file. Path may be local or global.
                May also be the text of a PDB file from the internet.
        """
        if '.pdb' in text: # means this is a filename
            self.pdbfile = text
            fname = text.split('.pdb')[0]
            f = open(fname + '.pdb', 'r')
            text = f.read()
            enter = '\n'
            f.close()
        else:
            enter = "\\n"
        # class attributes
        aas = set()
        hetatms = {}
        atoms = {}
        chains = {}
        missing_atoms = {}
        missing_aas = []
        conf = []
        bonds = {}
        # get R and Rfree values
        if "R VALUE            (WORKING SET)" in text:
            temp = text.split("R VALUE            (WORKING SET)")
            temp2 = temp[-1].split()
            if temp2[1] != 'NULL':
                R = float(temp2[1])
            else:
                R = -100
            if temp2[8] != 'NULL':
                Rfree = float(temp2[8])
            else:
                Rfree = 100
        else:
            temp = text.split("R VALUE          (WORKING SET, NO CUTOFF)")
            temp2 = temp[-1].split()
            if temp2[1] != 'NULL':
                R = float(temp2[1])
            else:
                R = -100
            if temp2[10] != 'NULL':
                Rfree = float(temp2[10])
            else:
                Rfree = 100
        temp = temp[1].split(enter)
        # start getting missing amino acids
        if "M RES C SSSEQI" in text:
            text = text.split("M RES C SSSEQI")
            want = text[-1]
            text = text[0].split(enter)
            split = text[-1]
            want = want.split(split)
            for line in want:
                if line == want[-1]:
                    text = line
                    line = line.split(enter)
                    line = line[0]
                    text = text.replace(line, '')
                l = line.split()
                if len(l) > 2:
                    a = AA3D(l[0], l[1], l[2])
                    missing_aas.append(a)
        # start getting missing atoms
        if "M RES CSSEQI  ATOMS" in text:
            text = text.split("M RES CSSEQI  ATOMS")
            want = text[-1]
            text = text[0].split(enter)
            split = text[-1]
            want = want.split(split)
            for line in want:
                if line == want[-1]: 
                    text = line
                    line = line.split(enter)
                    line = line[0]
                    text = text.replace(line, '')
                l = line.split()
                if len(l) > 2:
                    a = AA3D(l[0], l[1], l[2])
                    missing_atoms[a] = []
                    for atom in l[3:]:
                        if atom != enter and atom[0] in ['C', 'N', 'O', 'H']:
                            missing_atoms[a].append(atom3D(Sym=atom[0],
                                                           greek=atom))
        # start getting amino acids and heteroatoms
        if "ENDMDL" in text:
            text.split("ENDMDL")
            text = text[-2] + text[-1]
        text = text.split(enter)
        text = text[1:]
        for line in text:
            if line == text[-1]:
                text = line
                line = line.split(enter)
                line = line[0]
                text = text.replace(line, '')
            l = line.split()
            l_type = l[0]
            l = l[1:]
            if "ATOM" in l_type: # we are in an amino acid
                if len(l_type) > 4: # fixes buggy splitting at start
                    l = [l_type[4:]] + l
                l = pStripSpaces(l)
                a = AA3D(l[2], l[3], l[4], float(l[8]))
                if l[3] not in chains.keys():
                    chains[l[3]] = [] # initialize key of chain dictionary
                if int(float(l[8])) != 1 and a not in conf:
                    conf.append(a)
                if a not in chains[l[3]] and a not in conf:
                    chains[l[3]].append(a)
                aas.add(a)
                atom = atom3D(Sym=l[10], xyz=[l[5], l[6], l[7]], Tfactor=l[9],
                              occup=float(l[8]), greek=l[1])
                a.addAtom(atom, int(l[0])) # terminal Os may be missing
                atoms[int(l[0])] = atom
                a.setBonds()
                bonds.update(a.bonds)
                if a.prev != None:
                    bonds[a.n].add(a.prev.c)
                if a.next != None:
                    bonds[a.c].add(a.next.n)
            elif "HETATM" in l_type: # this is a heteroatom
                if len(l_type) > 6: # fixes buggy splitting at start
                    l = [l_type[6:]] + l
                l = pStripSpaces(l)
                hetatm = atom3D(Sym=l[-1], xyz = [l[5], l[6], l[7]], Tfactor=l[9],
                                occup=float(l[8]), greek=l[1])
                if (int(l[0]), hetatm) not in hetatms.keys():
                    hetatms[(int(l[0]), hetatm)] = [l[2], l[3]] # [cmpd name, chain]
                atoms[int(l[0])] = hetatm
            elif "CONECT" in l_type: # get extra connections
                if len(l_type) > 6: # fixes buggy splitting
                    l = [l_type[6:]] + l
                l2 = []
                wrong_turn = False
                for i in range(len(l)):
                    x = l[i]
                    wrong_turn = False
                    while x != '' and int(x) not in atoms.keys():
                        if int(x[:5]) in atoms.keys():
                            l2.append(x[:5])
                            x = x[5:]
                            wrong_turn = False
                        elif int(x[:4]) in atoms.keys():
                            l2.append(x[:4])
                            x = x[4:]
                            wrong_turn = False
                        elif (l2 == [] or wrong_turn) and int(x[:3]) in atoms.keys():
                            l2.append(x[:3])
                            x = x[3:]
                            wrong_turn = False
                        else: # made a wrong turn
                            wrong_turn = True
                            y = l2.pop()
                            l2.append(y[:-1])
                            x = y[-1] + x
                    if x != '':
                        l2.append(x)
                l = l2
                if l != [] and atoms[int(l[0])] not in bonds.keys():
                    bonds[atoms[int(l[0])]] = set()
                for i in l[1:]:
                    if i != '':
                        bonds[atoms[int(l[0])]].add(atoms[int(i)])
        # deal with conformations
        for i in range(len(conf)-1):
            if conf[i].chain == conf[i+1].chain and conf[i].id == conf[i+1].id:
                if conf[i].occup >= conf[i+1].occup:
                    chains[conf[i].chain].append(conf[i])
                    # pick chain with higher occupancy or the A chain if tie
                else:
                    chains[conf[i+1].chain].append(conf[i+1])
        self.setChains(chains)
        self.setAAs(aas)
        self.setAtoms(atoms)
        self.setHetatms(hetatms)
        self.setMissingAtoms(missing_atoms)
        self.setMissingAAs(missing_aas)
        self.setConf(conf)
        self.setR(R)
        self.setRfree(Rfree)
        self.setBonds(bonds)

    def fetch_pdb(self, pdbCode):
        """ API query to fetch a pdb and write it as a protein3D class instance

        Parameters
        ----------
            pdbCode : str
                code for protein, e.g. 1os7
        """
        remoteCode = pdbCode.upper()
        try:
            data = urllib.urlopen(
                'https://files.rcsb.org/view/' + remoteCode +
                '.pdb').read()
        except:
            print("warning: %s not found.\n"%pdbCode)
        else:
            try:
                self.readfrompdb(str(data))
                print("fetched: %s"%(pdbCode))
            except IOError:
                print('aborted')
            else:
                if len(data) == 0:
                    print("warning: %s not valid.\n"%pdbCode)

    def setBonds(self, bonds):
        """Sets the bonded atoms in the protein.
        This is effectively the molecular graph.
         
        Parameters
        ----------
            bonds : dictionary
                Keyed by atom3D atoms in the protein
                Valued by a set consisting of bonded atoms
        """
        self.bonds = bonds

    def readMetaData(self, pdbCode):
        """ API query to fetch XML data from a pdb and add its useful attributes
        to a protein3D class.
        
        Parameters
        ----------
            pdbCode : str
                code for protein, e.g. 1os7
        """
        try:
            start = 'https://files.rcsb.org/pub/pdb/validation_reports/' + pdbCode[1] + pdbCode[2]
            link = start + '/' + pdbCode + '/' + pdbCode + '_validation.xml'
            xml_doc = requests.get(link)
        except:
            print("warning: %s not found.\n"%pdbCode)
        else:
            try:
                ### We then use beautiful soup to read the XML doc. LXML is an XML reader. The soup object is what we then use to parse!
                soup = BeautifulSoup(xml_doc.content,'lxml-xml')
                
                ### We can then use methods of the soup object to find "tags" within the XML file. This is how we would extract sections. 
                ### This is an example of getting everything with a "sec" tag.
                body = soup.find_all('wwPDB-validation-information')
                entry = body[0].find_all("Entry")
                #print("got metadata: %s"%(pdbCode))
                if "DataCompleteness" not in entry[0].attrs.keys():
                    self.setDataCompleteness(0)
                    print("warning: %s has no DataCompleteness."%pdbCode)
                else:
                    self.setDataCompleteness(float(entry[0].attrs["DataCompleteness"]))
                if "percent-RSRZ-outliers" not in entry[0].attrs.keys():
                    self.setRSRZ(100)
                    print("warning: %s has no RSRZ.\n"%pdbCode)
                else:
                    self.setRSRZ(float(entry[0].attrs["percent-RSRZ-outliers"]))
                if "TwinL" not in entry[0].attrs.keys():
                    print("warning: %s has no TwinL."%pdbCode)
                    self.setTwinL(0)
                else:
                    self.setTwinL(float(entry[0].attrs["TwinL"]))
                if "TwinL2" not in entry[0].attrs.keys():
                    print("warning: %s has no TwinL2."%pdbCode)
                    self.setTwinL2(0)
                else:
                    self.setTwinL2(float(entry[0].attrs["TwinL2"]))
            except IOError:
                print('aborted')
            else:
                if xml_doc == None:
                    print("warning: %s not valid.\n"%pdbCode)

    def setDataCompleteness(self, DataCompleteness):
        """ Set DataCompleteness value of protein3D class.

        Parameters
        ----------
            DataCompleteness : float
                The desired new R value.
        """
        self.DataCompleteness = DataCompleteness

    def setTwinL(self, TwinL):
        """ Set TwinL score of protein3D class.

        Parameters
        ----------
            TwinL : float
                The desired new TwinL score.
        """
        self.TwinL = TwinL

    def setTwinL2(self, TwinL2):
        """ Set TwinL squared score of protein3D class.

        Parameters
        ----------
            TwinL2 : float
                The desired new TwinL squared score.
        """
        self.TwinL2 = TwinL2

