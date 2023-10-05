def get_limits_r():
    # Numbers copied from "../files/UpperLimit_ratio_log.C" made by Fanqiang. To be changed whenever the limits change (if ever).
    y_95 = [0.005, 0.005271532, 0.005557811, 0.005859636, 0.006177852, 0.006513349, 0.006867066, 0.007239992, 0.007633171, 0.008047701, 0.008484743, 0.00894552, 0.009431319, 0.009943501, 0.0104835, 0.01105282, 0.01165306, 0.0122859, 0.0129531, 0.01365654, 0.01439817, 0.01518009, 0.01600446, 0.01687361, 0.01778996, 0.01875607, 0.01977464, 0.02084853, 0.02198074, 0.02317444, 0.02443296, 0.02575983, 0.02715876, 0.02863365, 0.03018864, 0.03182808, 0.03355655, 0.03537889, 0.03730019, 0.03932584, 0.04146148, 0.04371311, 0.04608701, 0.04858984, 0.05122858, 0.05401062, 0.05694375, 0.06003616, 0.06329651, 0.06673392, 0.07035801, 0.0741789, 0.0782073, 0.08245446, 0.08693227, 0.09165325, 0.09663062, 0.1018783, 0.1074109, 0.113244, 0.1193939, 0.1258778, 0.1327138, 0.139921, 0.1475196, 0.1555309, 0.1639772, 0.1728822, 0.1822708, 0.1921693, 0.2026054, 0.2136081, 0.2252084, 0.2374387, 0.2503332, 0.2639279, 0.2782609, 0.2933723, 0.3093043, 0.3261015, 0.3438109, 0.3624821, 0.3821672, 0.4029213, 0.4248026, 0.4478721, 0.4721944, 0.4978377, 0.5248735, 0.5533775, 0.5834295, 0.6151134, 0.6485181, 0.6837368, 0.7208681, 0.7600159, 0.8012897, 0.8448049, 0.8906833, 0.9390532, 0.9900498, 0.9900498, 0.9390532, 0.8906833, 0.8448049, 0.8012897, 0.7600159, 0.7208681, 0.6837368, 0.6485181, 0.6151134, 0.5834295, 0.5533775, 0.5248735, 0.4978377, 0.4721944, 0.4478721, 0.4248026, 0.4029213, 0.3821672, 0.3624821, 0.3438109, 0.3261015, 0.3093043, 0.2933723, 0.2782609, 0.2639279, 0.2503332, 0.2374387, 0.2252084, 0.2136081, 0.2026054, 0.1921693, 0.1822708, 0.1728822, 0.1639772, 0.1555309, 0.1475196, 0.139921, 0.1327138, 0.1258778, 0.1193939, 0.113244, 0.1074109, 0.1018783, 0.09663062, 0.09165325, 0.08693227, 0.08245446, 0.0782073, 0.0741789, 0.07035801, 0.06673392, 0.06329651, 0.06003616, 0.05694375, 0.05401062, 0.05122858, 0.04858984, 0.04608701, 0.04371311, 0.04146148, 0.03932584, 0.03730019, 0.03537889, 0.03355655, 0.03182808, 0.03018864, 0.02863365, 0.02715876, 0.02575983, 0.02443296, 0.02317444, 0.02198074, 0.02084853, 0.01977464, 0.01875607, 0.01778996, 0.01687361, 0.01600446, 0.01518009, 0.01439817, 0.01365654, 0.0129531, 0.0122859, 0.01165306, 0.01105282, 0.0104835, 0.009943501, 0.009431319, 0.00894552, 0.008484743, 0.008047701, 0.007633171, 0.007239992, 0.006867066, 0.006513349, 0.006177852, 0.005859636, 0.005557811, 0.005271532, 0.005]
    x_95 = [26.7274, 26.5588, 26.485, 26.377, 26.234, 26.0911, 26.0543, 25.912, 25.746, 25.5951, 25.4538, 25.3822, 25.2182, 25.0778, 24.9844, 24.8449, 24.7057, 24.5667, 24.4884, 24.328, 24.1903, 24.123, 23.955, 23.8187, 23.6827, 23.5472, 23.412, 23.2772, 23.1627, 22.9585, 22.8257, 22.6933, 22.5615, 22.3597, 22.2549, 22.151, 21.951, 21.8631, 21.6644, 21.4665, 21.3808, 21.1844, 21.0603, 20.9059, 20.7838, 20.6314, 20.4395, 20.2612, 20.1149, 19.8949, 19.7476, 19.6016, 19.4889, 19.2736, 19.073, 18.9328, 18.7217, 18.5241, 18.3888, 18.2237, 18.0195, 17.8275, 17.5951, 17.4062, 17.2517, 17.066, 16.8435, 16.6295, 16.4125, 16.2035, 15.9963, 15.715, 15.5148, 15.2836, 15.0891, 14.8273, 14.5694, 14.3157, 14.0668, 13.7896, 13.535, 13.2452, 12.9728, 12.6636, 12.3319, 12.0066, 11.6825, 11.335, 10.9481, 10.5529, 10.1315, 9.6958, 9.1984, 8.6632, 8.1326, 7.4833, 6.7259, 5.9077, 4.773, 3.2971, 1.1258, 0.25, 0.25, 0.5, 0.5, 0.75, 0.875, 1, 1.125, 0.6875, 1.125, 1.2188, 1.3125, 1.4062, 1.4531, 1.8047, 1.9141, 2.0234, 2.0781, 2.1875, 2.2969, 2.3516, 2.4609, 2.7246, 2.8125, 2.9004, 2.9883, 3.0762, 3.1641, 3.3301, 3.3398, 3.5117, 3.7188, 3.8125, 3.875, 3.9688, 4.0312, 4.125, 4.1875, 4.2812, 4.4795, 4.5762, 4.6406, 4.7051, 4.7695, 4.8662, 5.0801, 5.1465, 5.2461, 5.3125, 5.3789, 5.4453, 5.5928, 5.6602, 5.7444, 5.8286, 5.9814, 6.0498, 6.1182, 6.1865, 6.3442, 6.4136, 6.4829, 6.5522, 6.6216, 6.7852, 6.8555, 6.9082, 6.9609, 7.1289, 7.2002, 7.2715, 7.3428, 7.4141, 7.5007, 7.624, 7.6963, 7.7686, 7.8408, 7.9131, 8.0029, 8.0757, 8.2031, 8.2397, 8.313, 8.3862, 8.4595, 8.5327, 8.6836, 8.7578, 8.7949, 8.8691, 8.9434, 9.0176, 9.1143, 9.189, 9.2866, 9.3618, 9.437, 9.4746, 9.5498, 9.6497]
    (y_95_high, y_95_low) = y_95[:len(y_95)/2], y_95[len(y_95)/2:]
    (x_95_high, x_95_low) = x_95[:len(x_95)/2], x_95[len(x_95)/2:]
    y_68 = [0.005, 0.005271532, 0.005557811, 0.005859636, 0.006177852, 0.006513349, 0.006867066, 0.007239992, 0.007633171, 0.008047701, 0.008484743, 0.00894552, 0.009431319, 0.009943501, 0.0104835, 0.01105282, 0.01165306, 0.0122859, 0.0129531, 0.01365654, 0.01439817, 0.01518009, 0.01600446, 0.01687361, 0.01778996, 0.01875607, 0.01977464, 0.02084853, 0.02198074, 0.02317444, 0.02443296, 0.02575983, 0.02715876, 0.02863365, 0.03018864, 0.03182808, 0.03355655, 0.03537889, 0.03730019, 0.03932584, 0.04146148, 0.04371311, 0.04608701, 0.04858984, 0.05122858, 0.05401062, 0.05694375, 0.06003616, 0.06329651, 0.06673392, 0.07035801, 0.0741789, 0.0782073, 0.08245446, 0.08693227, 0.09165325, 0.09663062, 0.1018783, 0.1074109, 0.113244, 0.1193939, 0.1258778, 0.1327138, 0.139921, 0.1475196, 0.1555309, 0.1639772, 0.1728822, 0.1822708, 0.1921693, 0.2026054, 0.2136081, 0.2252084, 0.2374387, 0.2503332, 0.2639279, 0.2782609, 0.2933723, 0.3093043, 0.3261015, 0.3438109, 0.3624821, 0.3821672, 0.4029213, 0.4248026, 0.4478721, 0.4721944, 0.4978377, 0.5248735, 0.5533775, 0.5834295, 0.6151134, 0.6485181, 0.6837368, 0.7208681, 0.7600159, 0.8012897, 0.8448049, 0.8906833, 0.9390532, 0.9900498, 0.9900498, 0.9390532, 0.8906833, 0.8448049, 0.8012897, 0.7600159, 0.7208681, 0.6837368, 0.6485181, 0.6151134, 0.5834295, 0.5533775, 0.5248735, 0.4978377, 0.4721944, 0.4478721, 0.4248026, 0.4029213, 0.3821672, 0.3624821, 0.3438109, 0.3261015, 0.3093043, 0.2933723, 0.2782609, 0.2639279, 0.2503332, 0.2374387, 0.2252084, 0.2136081, 0.2026054, 0.1921693, 0.1822708, 0.1728822, 0.1639772, 0.1555309, 0.1475196, 0.139921, 0.1327138, 0.1258778, 0.1193939, 0.113244, 0.1074109, 0.1018783, 0.09663062, 0.09165325, 0.08693227, 0.08245446, 0.0782073, 0.0741789, 0.07035801, 0.06673392, 0.06329651, 0.06003616, 0.05694375, 0.05401062, 0.05122858, 0.04858984, 0.04608701, 0.04371311, 0.04146148, 0.03932584, 0.03730019, 0.03537889, 0.03355655, 0.03182808, 0.03018864, 0.02863365, 0.02715876, 0.02575983, 0.02443296, 0.02317444, 0.02198074, 0.02084853, 0.01977464, 0.01875607, 0.01778996, 0.01687361, 0.01600446, 0.01518009, 0.01439817, 0.01365654, 0.0129531, 0.0122859, 0.01165306, 0.01105282, 0.0104835, 0.009943501, 0.009431319, 0.00894552, 0.008484743, 0.008047701, 0.007633171, 0.007239992, 0.006867066, 0.006513349, 0.006177852, 0.005859636, 0.005557811, 0.005271532, 0.005]
    x_68 = [20.9055, 20.7649, 20.6314, 20.58, 20.4468, 20.3136, 20.2628, 20.1299, 19.9907, 19.947, 19.8145, 19.682, 19.544, 19.412, 19.2743, 19.1427, 19.0113, 18.8799, 18.8371, 18.7011, 18.5703, 18.4396, 18.3931, 18.2629, 18.1329, 18.0031, 17.8734, 17.7439, 17.5667, 17.4379, 17.3094, 17.1811, 17.053, 16.9252, 16.7962, 16.6675, 16.5408, 16.3701, 16.244, 16.1183, 15.9487, 15.8237, 15.6991, 15.5306, 15.4069, 15.2393, 15.1166, 14.9505, 14.8288, 14.6632, 14.4982, 14.3337, 14.2141, 14.0509, 13.8877, 13.7258, 13.5647, 13.4027, 13.2431, 13.04, 12.882, 12.7228, 12.5224, 12.3644, 12.1661, 12.0096, 11.8137, 11.6147, 11.4213, 11.2245, 11.0285, 10.8398, 10.6394, 10.4102, 10.2197, 9.9861, 9.7542, 9.5243, 9.2964, 9.0261, 8.8349, 8.583, 8.3077, 8.0332, 7.7488, 7.4373, 7.1732, 6.8224, 6.5145, 6.1758, 5.8429, 5.424, 5.0124, 4.6664, 4.1865, 3.7194, 3.1762, 2.5964, 1.9095, 1.3174, 0.8132, 0.375, 0.375, 0.75, 0.75, 1.125, 1.0938, 1.25, 1.4062, 1.7188, 1.5938, 1.7266, 2.1328, 2.2852, 2.3613, 2.6748, 2.8369, 2.8364, 3.0801, 3.2422, 3.4043, 3.4854, 3.6475, 3.7861, 4.0078, 4.1331, 4.2583, 4.3835, 4.5088, 4.6594, 4.7593, 4.9136, 4.9971, 5.123, 5.3281, 5.457, 5.543, 5.6719, 5.7578, 5.8867, 6.0575, 6.1882, 6.2754, 6.3625, 6.5197, 6.6519, 6.761, 6.8494, 6.9819, 7.0703, 7.2328, 7.3221, 7.4608, 7.5507, 7.6631, 7.7754, 7.9174, 8.0079, 8.0984, 8.1889, 8.3338, 8.4249, 8.516, 8.607, 8.7396, 8.8466, 8.9383, 9.049, 9.1604, 9.227, 9.3193, 9.4543, 9.547, 9.6397, 9.7602, 9.8407, 9.934, 10.0273, 10.1653, 10.2589, 10.3377, 10.4317, 10.5128, 10.6052, 10.6995, 10.7938, 10.888, 10.9823, 11.0973, 11.1921, 11.2396, 11.3344, 11.4293, 11.5241, 11.6069, 11.702, 11.7854, 11.8808, 11.9763, 12.073, 12.1688, 12.2532]
    (y_68_high, y_68_low) = y_68[:len(y_68)/2], y_68[len(y_68)/2:]
    (x_68_high, x_68_low) = x_68[:len(x_68)/2], x_68[len(x_68)/2:]
    y_exp = [0.005, 0.005271532, 0.005557811, 0.005859636, 0.006177852, 0.006513349, 0.006867066, 0.007239992, 0.007633171, 0.008047701, 0.008484743, 0.00894552, 0.009431319, 0.009943501, 0.0104835, 0.01105282, 0.01165306, 0.0122859, 0.0129531, 0.01365654, 0.01439817, 0.01518009, 0.01600446, 0.01687361, 0.01778996, 0.01875607, 0.01977464, 0.02084853, 0.02198074, 0.02317444, 0.02443296, 0.02575983, 0.02715876, 0.02863365, 0.03018864, 0.03182808, 0.03355655, 0.03537889, 0.03730019, 0.03932584, 0.04146148, 0.04371311, 0.04608701, 0.04858984, 0.05122858, 0.05401062, 0.05694375, 0.06003616, 0.06329651, 0.06673392, 0.07035801, 0.0741789, 0.0782073, 0.08245446, 0.08693227, 0.09165325, 0.09663062, 0.1018783, 0.1074109, 0.113244, 0.1193939, 0.1258778, 0.1327138, 0.139921, 0.1475196, 0.1555309, 0.1639772, 0.1728822, 0.1822708, 0.1921693, 0.2026054, 0.2136081, 0.2252084, 0.2374387, 0.2503332, 0.2639279, 0.2782609, 0.2933723, 0.3093043, 0.3261015, 0.3438109, 0.3624821, 0.3821672, 0.4029213, 0.4248026, 0.4478721, 0.4721944, 0.4978377, 0.5248735, 0.5533775, 0.5834295, 0.6151134, 0.6485181, 0.6837368, 0.7208681, 0.7600159, 0.8012897, 0.8448049, 0.8906833, 0.9390532, 0.9900498]
    x_exp = [15.9375, 15.875, 15.75, 15.6875, 15.5625, 15.4375, 15.375, 15.25, 15.1875, 15.0625, 14.9375, 14.8125, 14.75, 14.625, 14.5625, 14.4375, 14.3125, 14.1875, 14.0625, 14, 13.875, 13.75, 13.6875, 13.5625, 13.4375, 13.3125, 13.1875, 13.0625, 13, 12.875, 12.75, 12.625, 12.5, 12.375, 12.2812, 12.1875, 12.0625, 11.9375, 11.8125, 11.6875, 11.5625, 11.4375, 11.3125, 11.1875, 11.0625, 10.9375, 10.8125, 10.6562, 10.5, 10.375, 10.25, 10.125, 10, 9.875, 9.6875, 9.5625, 9.4375, 9.25, 9.125, 9, 8.875, 8.6875, 8.5625, 8.375, 8.25, 8.0625, 7.9375, 7.75, 7.625, 7.4375, 7.25, 7.125, 6.875, 6.75, 6.5625, 6.375, 6.1875, 6, 5.8125, 5.625, 5.375, 5.25, 5, 4.75, 4.625, 4.375, 4.125, 3.875, 3.75, 3.5, 3.25, 3, 2.75, 2.25, 2, 1.75, 1.5, 1, 1, 0.5, 0.5]
    x_obs = [10.9389, 10.8679, 10.7812, 10.6943, 10.6071, 10.5196, 10.4318, 10.3438, 10.2305, 10.1447, 10.078, 9.9677, 9.8954, 9.7986, 9.7113, 9.6238, 9.506, 9.4148, 9.3234, 9.2316, 9.1395, 9.0471, 8.9548, 8.8772, 8.7838, 8.6901, 8.596, 8.5016, 8.4068, 8.3117, 8.1924, 8.0996, 7.9963, 7.9182, 7.8187, 7.7242, 7.6136, 7.5205, 7.4217, 7.3225, 7.2229, 7.1228, 7.0222, 6.9207, 6.804, 6.7021, 6.5997, 6.4969, 6.3935, 6.2897, 6.1831, 6.0776, 5.972, 5.8655, 5.7584, 5.6566, 5.5483, 5.4402, 5.3315, 5.2223, 5.1124, 5.002, 4.8751, 4.7634, 4.6511, 4.5382, 4.4245, 4.3103, 4.2, 4.0757, 3.9665, 3.8498, 3.7237, 3.6044, 3.4904, 3.3634, 3.2417, 3.1192, 2.9959, 2.8749, 2.7499, 2.6239, 2.4971, 2.3693, 2.2406, 2.111, 1.9798, 1.8481, 1.7154, 1.5793, 1.4443, 1.3082, 1.171, 1.0333, 0.8938, 0.7518, 0.6098, 0.4665, 0.3219, 0.1755, 0.0278]
    y_obs = [0.005, 0.005271532, 0.005557811, 0.005859636, 0.006177852, 0.006513349, 0.006867066, 0.007239992, 0.007633171, 0.008047701, 0.008484743, 0.00894552, 0.009431319, 0.009943501, 0.0104835, 0.01105282, 0.01165306, 0.0122859, 0.0129531, 0.01365654, 0.01439817, 0.01518009, 0.01600446, 0.01687361, 0.01778996, 0.01875607, 0.01977464, 0.02084853, 0.02198074, 0.02317444, 0.02443296, 0.02575983, 0.02715876, 0.02863365, 0.03018864, 0.03182808, 0.03355655, 0.03537889, 0.03730019, 0.03932584, 0.04146148, 0.04371311, 0.04608701, 0.04858984, 0.05122858, 0.05401062, 0.05694375, 0.06003616, 0.06329651, 0.06673392, 0.07035801, 0.0741789, 0.0782073, 0.08245446, 0.08693227, 0.09165325, 0.09663062, 0.1018783, 0.1074109, 0.113244, 0.1193939, 0.1258778, 0.1327138, 0.139921, 0.1475196, 0.1555309, 0.1639772, 0.1728822, 0.1822708, 0.1921693, 0.2026054, 0.2136081, 0.2252084, 0.2374387, 0.2503332, 0.2639279, 0.2782609, 0.2933723, 0.3093043, 0.3261015, 0.3438109, 0.3624821, 0.3821672, 0.4029213, 0.4248026, 0.4478721, 0.4721944, 0.4978377, 0.5248735, 0.5533775, 0.5834295, 0.6151134, 0.6485181, 0.6837368, 0.7208681, 0.7600159, 0.8012897, 0.8448049, 0.8906833, 0.9390532, 0.9900498]

    return (zip(x_obs, y_obs), zip(x_exp, y_exp), zip(x_68_low, y_68_low), zip(x_68_high, y_68_high), zip(x_95_low, y_95_low), zip(x_95_high, y_95_high))