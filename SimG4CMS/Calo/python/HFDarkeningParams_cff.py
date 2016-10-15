import FWCore.ParameterSet.Config as cms

# These parameters are for modeling HF fiber radiation damage and resulting darkening.
#

HFDarkeningParameterBlock = cms.PSet(
  doseLayerDepth = cms.vdouble(
   *[
      #Depth 0: 1115. - 1120. cm
      1089.82672, 430.74392, 316.10688, 198.15032, 133.66563, 92.69111, 63.84813, 45.19114, 30.58487, 21.00216, 12.95498, 5.45337, 0.4561,

      #Depth 1: 1120. - 1125. cm
      1662.5198, 720.46188, 590.62325, 429.82339, 282.39827, 186.10415, 120.91404, 79.7575, 51.45268, 32.77064, 17.95221, 6.07963, 0.32217,

      #Depth 2: 1125. - 1130. cm
      1511.32547, 708.54584, 605.99167, 495.66233, 322.7595, 207.45045, 129.52382, 81.14806, 51.51518, 30.54072, 15.88012, 5.04118, 0.29664,

      #Depth 3: 1130. - 1135. cm
      1116.32013, 574.64553, 480.43324, 409.50734, 266.19544, 168.75067, 104.73339, 63.34151, 39.86312, 22.97123, 11.55581, 3.68691, 0.29269,

      #Depth 4: 1135. - 1140. cm
      800.2853, 446.24557, 359.63381, 299.49924, 195.82939, 124.59848, 77.29738, 46.00855, 28.52664, 16.43682, 8.18736, 2.67659, 0.28495,

      #Depth 5: 1140. - 1145. cm
      612.20865, 355.40823, 276.97155, 219.61391, 144.74911, 93.28837, 58.51092, 35.08033, 21.47604, 12.29605, 6.16789, 2.12159, 0.27661,

      #Depth 6: 1145. - 1150. cm
      504.87242, 295.98737, 226.62047, 171.60095, 113.64104, 73.25579, 46.67008, 28.26167, 17.42322, 9.59497, 5.02061, 1.76107, 0.28319,

      #Depth 7: 1150. - 1155. cm
      439.39199, 252.97218, 191.68489, 143.625, 94.50761, 61.46441, 38.8409, 23.68308, 14.69344, 7.98846, 4.18968, 1.52224, 0.28134,

      #Depth 8: 1155. - 1160. cm
      400.29378, 223.24941, 165.77051, 124.7812, 80.69381, 52.80131, 33.61851, 20.37541, 12.53798, 7.00165, 3.62975, 1.38213, 0.26854,

      #Depth 9: 1160. - 1165. cm
      367.42369, 201.03927, 146.05054, 108.58533, 70.57479, 46.07769, 29.50341, 17.87488, 10.91036, 6.12909, 3.21967, 1.25779, 0.25991,

      #Depth 10: 1165. - 1170. cm
      329.90783, 181.01783, 128.96837, 94.23766, 62.25865, 40.17712, 25.787, 16.16342, 9.46068, 5.41504, 2.81787, 1.13574, 0.2552,

      #Depth 11: 1170. - 1175. cm
      297.89937, 163.50102, 116.23316, 83.18498, 55.46801, 35.4402, 22.19535, 14.31524, 8.21573, 4.69969, 2.50891, 1.02412, 0.24167,

      #Depth 12: 1175. - 1180. cm
      272.26937, 147.64059, 104.98641, 72.05083, 48.82044, 31.0403, 19.50233, 12.6764, 7.06137, 4.11597, 2.14441, 0.93751, 0.22596,

      #Depth 13: 1180. - 1185. cm
      249.14766, 131.41739, 93.58642, 63.48269, 42.64957, 27.27858, 17.06007, 10.93655, 6.0742, 3.64389, 1.88608, 0.84005, 0.22152,

      #Depth 14: 1185. - 1190. cm
      226.99648, 116.79805, 82.04924, 55.57321, 37.60804, 23.8003, 15.12373, 9.59352, 5.28989, 3.16053, 1.65219, 0.74207, 0.20077, 

      #Depth 15: 1190. - 1195. cm
      207.71168, 102.36457, 72.48674, 48.65921, 32.76174, 20.67312, 13.34669, 8.26618, 4.70605, 2.83691, 1.4537, 0.66819, 0.18359, 

      #Depth 16: 1195. - 1200. cm
      190.76744, 90.63229, 63.78329, 42.44254, 29.13414, 18.03128, 11.85671, 7.31436, 4.02773, 2.45523, 1.26776, 0.58725, 0.18424, 

      #Depth 17: 1200. - 1205. cm
      178.21925, 81.88379, 57.03839, 37.22998, 24.90575, 16.00853, 10.46851, 6.57621, 3.50808, 2.16873, 1.12083, 0.51166, 0.16961, 

      #Depth 18: 1205. - 1210. cm
      169.78037, 73.27996, 50.41326, 32.70998, 21.81181, 14.35241, 9.10403, 5.68453, 3.12387, 1.85718, 1.00079, 0.44896, 0.15281, 

      #Depth 19: 1210. - 1215. cm
      162.59574, 66.56493, 43.80282, 28.71882, 18.54598, 12.41703, 7.92044, 4.83008, 2.7721, 1.5999, 0.86314, 0.37843, 0.13917, 

      #Depth 20: 1215. - 1220. cm
      158.18666, 59.65062, 38.26824, 24.43072, 16.41785, 10.76117, 6.82664, 4.16347, 2.44142, 1.4081, 0.76763, 0.34437, 0.13022, 

      #Depth 21: 1220. - 1225. cm
      150.13577, 53.33958, 34.18336, 21.62898, 14.27769, 9.37977, 5.89517, 3.6894, 2.16963, 1.18001, 0.6684, 0.3026, 0.12473, 

      #Depth 22: 1225. - 1230. cm
      144.67223, 50.31662, 31.59758, 19.6959, 12.52437, 8.238, 5.24535, 3.24359, 1.87348, 1.03871, 0.58652, 0.26342, 0.11119, 

      #Depth 23: 1230. - 1235. cm
      136.98616, 45.60196, 28.92017, 17.66014, 10.97254, 7.36928, 4.53706, 2.80986, 1.65139, 0.89451, 0.49702, 0.22609, 0.09944, 

      #Depth 24: 1235. - 1240. cm
      132.02552, 41.30027, 26.02537, 15.64933, 9.94368, 6.42028, 3.90554, 2.47866, 1.47008, 0.8163, 0.4273, 0.20726, 0.08913, 

      #Depth 25: 1240. - 1245. cm
      129.36562, 38.47342, 22.91677, 13.76957, 8.7923, 5.50856, 3.44091, 2.1008, 1.26223, 0.69379, 0.38069, 0.17625, 0.07827, 

      #Depth 26: 1245. - 1250. cm
      127.78974, 35.60432, 20.49112, 12.18895, 7.72131, 4.74755, 2.91111, 1.82576, 1.12159, 0.5982, 0.34274, 0.15086, 0.06907, 

      #Depth 27: 1250. - 1255. cm
      124.07525, 33.24689, 18.53004, 10.8278, 6.97957, 4.0935, 2.46634, 1.63918, 0.97005, 0.50387, 0.29183, 0.13397, 0.05972, 

      #Depth 28: 1255. - 1260. cm
      122.43579, 31.75832, 17.29619, 9.83185, 6.05158, 3.55086, 2.21287, 1.39251, 0.82526, 0.45851, 0.24181, 0.11202, 0.05703, 

      #Depth 29: 1260. - 1265. cm
      121.61748, 30.44129, 16.17414, 8.93093, 5.44004, 3.08939, 1.905, 1.21964, 0.71349, 0.40864, 0.21309, 0.09492, 0.0535, 

      #Depth 30: 1265. - 1270. cm
      120.67237, 29.01728, 14.54675, 8.14094, 4.8913, 2.67395, 1.62092, 1.04203, 0.59052, 0.32811, 0.19279, 0.08357, 0.05398, 

      #Depth 31: 1270. - 1275. cm
      119.15188, 27.50844, 13.95703, 7.2369, 4.21802, 2.32933, 1.45072, 0.87999, 0.5075, 0.29475, 0.16882, 0.08048, 0.05269, 

      #Depth 32: 1275. - 1280. cm
      114.6991, 25.85671, 12.93505, 6.70319, 3.86928, 2.3702, 1.37139, 0.87629, 0.50304, 0.29778, 0.18413, 0.11134, 0.07111, 

      #Depth 33: 1280. - 1280. cm
      114.6991, 25.85671, 12.93505, 6.70319, 3.86928, 2.3702, 1.37139, 0.87629, 0.50304, 0.29778, 0.18413, 0.11134, 0.07111
     ])
)