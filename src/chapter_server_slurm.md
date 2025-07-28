# Setting Up a SLURM Cluster

This guide explains how to set up a SLURM workload manager on multiple servers, with one acting as the control/login node and the others as compute nodes.

## Prerequisites

- All servers run a compatible Linux distribution (e.g., Ubuntu, CentOS).
- All nodes have static IP addresses and can communicate over the network.
- Passwordless SSH is set up from the control node to all compute nodes.
- All nodes have synchronized clocks (e.g., via `chrony` or `ntpd`).

## 1. Install SLURM on All Nodes

On all nodes (control and compute):

```bash
sudo apt update
sudo apt install slurm-wlm munge
```
Or, for CentOS/RHEL:
```bash
sudo yum install slurm munge
```

## 2. Configure MUNGE Authentication

- On the control node, generate a MUNGE key:
    ```bash
    sudo /usr/sbin/create-munge-key
    ```
- Copy `/etc/munge/munge.key` to the same location on all compute nodes.
- Set correct permissions:
    ```bash
    sudo chown munge:munge /etc/munge/munge.key
    sudo chmod 400 /etc/munge/munge.key
    ```

## 3. Configure SLURM

- On the control node, generate a basic `slurm.conf` using [SLURM Configurator](https://slurm.schedmd.com/configurator.html).
- Example minimal `slurm.conf`:
    ```ini
    ControlMachine=control-node-hostname
    NodeName=compute[1-4] CPUs=4 State=UNKNOWN
    PartitionName=debug Nodes=ALL Default=YES MaxTime=INFINITE State=UP
    ```
- Copy `slurm.conf` to `/etc/slurm-llnl/slurm.conf` (or `/etc/slurm/slurm.conf`) on all nodes.

### Setting Up Multiple Partitions

Partitions in SLURM allow you to group nodes with different characteristics and apply different policies. Here are common partition configurations:

#### Example Configuration with Multiple Partitions:

```ini
# Control node configuration
ControlMachine=control-node-hostname
ControlAddr=192.168.1.10

# Node definitions
NodeName=compute[1-2] CPUs=8 RealMemory=16000 Gres=gpu:1 State=UNKNOWN
NodeName=compute[3-6] CPUs=4 RealMemory=8000 State=UNKNOWN
NodeName=compute[7-8] CPUs=16 RealMemory=32000 State=UNKNOWN

# Partition definitions
PartitionName=gpu Nodes=compute[1-2] Default=NO MaxTime=48:00:00 State=UP
PartitionName=normal Nodes=compute[3-6] Default=YES MaxTime=24:00:00 State=UP
PartitionName=bigmem Nodes=compute[7-8] Default=NO MaxTime=72:00:00 State=UP
PartitionName=all Nodes=compute[1-8] Default=NO MaxTime=INFINITE State=UP
```

#### Partition Configuration Options:

- **PartitionName**: Unique name for the partition
- **Nodes**: Which nodes belong to this partition (can overlap between partitions)
- **Default**: Whether this is the default partition (only one should be YES)
- **MaxTime**: Maximum job runtime (format: days-hours:minutes:seconds)
- **State**: Partition state (UP, DOWN, DRAIN, INACTIVE)
- **AllowGroups**: Restrict access to specific user groups
- **DenyGroups**: Deny access to specific user groups
- **MaxNodes**: Maximum nodes a job can use in this partition
- **MinNodes**: Minimum nodes required for jobs in this partition
- **Priority**: Partition priority (higher = more priority)

#### Advanced Partition Examples:

```ini
# High-priority partition for urgent jobs
PartitionName=urgent Nodes=compute[1-4] Default=NO MaxTime=4:00:00 State=UP Priority=1000

# Long-running jobs partition
PartitionName=long Nodes=compute[5-8] Default=NO MaxTime=168:00:00 State=UP Priority=1

# Interactive partition with short time limits
PartitionName=interactive Nodes=compute[1-2] Default=NO MaxTime=2:00:00 State=UP

# Research group specific partition
PartitionName=biogroup Nodes=compute[3-6] Default=NO AllowGroups=biolab MaxTime=48:00:00 State=UP
```

#### Using Partitions:

Users can specify partitions when submitting jobs:

```bash
# Submit to specific partition
sbatch -p gpu script.sh
sbatch --partition=bigmem job.sh

# Submit to multiple partitions (job runs on first available)
sbatch -p gpu,normal script.sh

# Check available partitions
sinfo
sinfo -p gpu
```

#### Best Practices for Partitions:

1. **Create purpose-specific partitions**: GPU nodes, high-memory nodes, etc.
2. **Set appropriate time limits**: Shorter limits for interactive work, longer for batch jobs
3. **Use priorities**: Higher priority for urgent/paying users
4. **Group restrictions**: Limit expensive resources to specific research groups
5. **Overlapping partitions**: Nodes can belong to multiple partitions for flexibility

## 4. Start and Enable Services

On all nodes:
```bash
sudo systemctl enable munge
sudo systemctl start munge
```

On the control node:
```bash
sudo systemctl enable slurmctld
sudo systemctl start slurmctld
```

On compute nodes:
```bash
sudo systemctl enable slurmd
sudo systemctl start slurmd
```

## 5. Verify the Cluster

- On the control node, check node status:
    ```bash
    sinfo
    scontrol show nodes
    ```
- Submit a test job:
    ```bash
    sbatch --wrap="hostname"
    ```

## 6. (Optional) Set Up SLURM Accounting

- Install and configure `slurmdbd` and a database (e.g., MariaDB) if accounting is needed.

## References

- [SLURM Quick Start Guide](https://slurm.schedmd.com/quickstart.html)
- [SLURM Documentation](https://slurm.schedmd.com/documentation.html)
