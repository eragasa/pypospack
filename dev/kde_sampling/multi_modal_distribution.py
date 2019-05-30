import os,shutil
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

plot_directory = 'kde'
plot_format_type = 'svg'
plot_dpi = 1200
multi_modal_pdf_name = os.path.join(plot_directory,'multi_modal_histogram_pdf.svg')
multi_modal_cdf_name = os.path.join(plot_directory,'multi_modal_histogram_cdf.svg')

n_samples = 1000
n_bins = 30

# create plot directory
if os.path.isdir(plot_directory):
    shutil.rmtree(plot_directory)
os.mkdir(plot_directory)

mu_1 = -5.0
sigma_1 = 2

mu_2 = 0.0
sigma_2 = 1.0

mu_3 = 5.0
sigma_3 = 2

norm_rv_1 = norm(loc=mu_1,scale=sigma_1)
norm_rv_2 = norm(loc=mu_2,scale=sigma_2)
norm_rv_3 = norm(loc=mu_3,scale=sigma_3)

# determin x_limits
x_min = min(norm_rv_1.ppf(0.0001),
            norm_rv_2.ppf(0.0001),
            norm_rv_3.ppf(0.0001))
x_max = max(norm_rv_1.ppf(0.9999),
            norm_rv_2.ppf(0.9999),
            norm_rv_3.ppf(0.9999))

x_min = -10
x_max = 10

x = np.linspace(x_min,x_max,500)

import scipy.integrate as integrate

pdf_1 = norm_rv_1.pdf(x)
pdf_2 = norm_rv_2.pdf(x)
pdf_3 = norm_rv_3.pdf(x)

c_pdf_1 = 1
c_pdf_2 = 2
c_pdf_3 = 1

pdf = c_pdf_1*pdf_1 \
      + c_pdf_2*pdf_2 \
      + c_pdf_3*pdf_3
pdf_normalization_factor = integrate.quad(
        lambda x: c_pdf_1*norm_rv_1.pdf(x) \
                  + c_pdf_2*norm_rv_2.pdf(x) \
                  + c_pdf_3*norm_rv_3.pdf(x),
        x_min,
        x_max)
pdf = pdf/pdf_normalization_factor[0]

#calculating the CDF by using the trapazoidal rule
print('calculating the CDF')
print('len(x)={}'.format(len(x)))
cdf = [np.trapz(pdf[:i],x[:i]) for i in range(len(x))]
cdf = np.array(cdf)

base_n_samples = int(n_samples/(c_pdf_1+c_pdf_2+c_pdf_3))
sample_rv = np.concatenate(
        (
            norm_rv_1.rvs(size=c_pdf_1*base_n_samples),
            norm_rv_2.rvs(size=c_pdf_2*base_n_samples),
            norm_rv_3.rvs(size=c_pdf_3*base_n_samples)
        )
)

fig, ax = plt.subplots(1,1)
ax.hist(sample_rv,
        bins=n_bins,
        density=True,
        histtype='stepfilled',
        alpha=0.2)
ax.plot(x,pdf)
ax.set_xlabel('x')
ax.set_ylabel('density')
ax.set_xlim(left=x_min,right=x_max)

fig.savefig(
        multi_modal_pdf_name,
        format=plot_format_type,
        dpi=plot_dpi)
plt.close('all')

# create cdf plot
print('creating cdf plot')
fig, ax = plt.subplots(1,1)


ax.hist(sample_rv,
        bins=n_bins,
        density=True,
        histtype='step',
        cumulative=True,
        color='black',
        label='sampling')


ax.plot(x,cdf)
ax.set_xlabel('x')
ax.set_ylabel('density')
ax.set_xlim(left=x_min,right=x_max)
fig.savefig(
        multi_modal_cdf_name,
        format=plot_format_type,
        dpi=plot_dpi)
plt.close(fig)

