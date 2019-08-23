
import matplotlib.pyplot as plt
import sys
sys.path.append('C:/Users/.../RMC_plotter/')
import RMC_plotter as rp


path = r"C:\RMCProfile\RMCProfile_package\tutorial\ex_6\rmc_neutron\run"
stem_name = "gapo4_neutron"
pdf_type = "gr"
recip_type = "sq"

fig, c, m = rp.rmcplot(path, stem_name, pdf_type, recip_type)

plt.show()




