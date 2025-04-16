import random

# Experimental migration rates
exp_migration_rates = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10])

# Compute theoretical F_ST values and add some noise
exp_F_ST_values = 1 / (1 + 4 * exp_migration_rates)
exp_F_ST_values = [max(0, min(1, fst + random.uniform(-0.05, 0.05))) for fst in exp_F_ST_values]  # Adding small noise

# Plot
plt.figure(figsize=(7, 5))
plt.plot(Nm_values, F_ST_values, label=r"$F_{ST} = \frac{1}{1 + 4Nm}$", color='b', lw=2)
plt.scatter(exp_migration_rates, exp_F_ST_values, color='r', label="Experimental Data", zorder=3)
plt.xscale("log")  # Log scale for Nm
plt.ylim(0, 1)
plt.xlabel("Migration Rate (Nm)", fontsize=12)
plt.ylabel(r"$F_{ST}$", fontsize=12)
plt.title("Relationship Between Migration Rate and F_ST", fontsize=14)
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.legend()
plt.show()
