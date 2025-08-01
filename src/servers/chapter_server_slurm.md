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
    ClusterName=control-node-hostname
    ControlMachine=control-node-hostname
    NodeName=compute[1-4] CPUs=4 State=UNKNOWN
    PartitionName=debug Nodes=ALL Default=YES MaxTime=INFINITE State=UP
    ```

**Important Configuration Notes:**
- `ClusterName`: Logical name of your SLURM cluster (can be descriptive)
- `ControlMachine`: Hostname of the control node (must match actual hostname)
- For single control node setups, these are typically the same
- Example: If your control node hostname is `mnemosyne2`, both should be:
    ```ini
    ClusterName=mnemosyne2
    ControlMachine=mnemosyne2
    ```

- Copy `slurm.conf` to `/etc/slurm-llnl/slurm.conf` (or `/etc/slurm/slurm.conf`) on all nodes.

### Setting Up Multiple Partitions

Partitions in SLURM allow you to group nodes with different characteristics and apply different policies. Here are common partition configurations:

#### Example Configuration with Multiple Partitions:

```ini
# Control node configuration
ControlMachine=control-node-hostname
ControlAddr=192.168.1.10

# Node definitions with sequential naming
NodeName=compute[1-2] CPUs=8 RealMemory=16000 Gres=gpu:1 State=UNKNOWN
NodeName=compute[3-6] CPUs=4 RealMemory=8000 State=UNKNOWN
NodeName=compute[7-8] CPUs=16 RealMemory=32000 State=UNKNOWN

# Partition definitions
PartitionName=gpu Nodes=compute[1-2] Default=NO MaxTime=48:00:00 State=UP
PartitionName=normal Nodes=compute[3-6] Default=YES MaxTime=24:00:00 State=UP
PartitionName=bigmem Nodes=compute[7-8] Default=NO MaxTime=72:00:00 State=UP
PartitionName=all Nodes=compute[1-8] Default=NO MaxTime=INFINITE State=UP
```

#### Example Configuration with Mixed/Arbitrary Node Names:

SLURM supports any node naming convention. Here's how to configure nodes with different naming patterns:

```ini
# Control node configuration
ControlMachine=control-node-hostname
ControlAddr=192.168.1.10

# Individual node definitions (different naming patterns)
NodeName=c1 CPUs=8 RealMemory=16000 Gres=gpu:1 State=UNKNOWN
NodeName=cc1 CPUs=4 RealMemory=8000 State=UNKNOWN
NodeName=cc2 CPUs=4 RealMemory=8000 State=UNKNOWN
NodeName=cr5 CPUs=16 RealMemory=32000 State=UNKNOWN
NodeName=txt455 CPUs=8 RealMemory=12000 State=UNKNOWN

# You can also group nodes with similar specs
NodeName=worker[1-3] CPUs=8 RealMemory=16000 State=UNKNOWN
NodeName=gpu-node1,gpu-node2 CPUs=12 RealMemory=24000 Gres=gpu:2 State=UNKNOWN

# Partition definitions with mixed node names
PartitionName=gpu Nodes=c1,gpu-node1,gpu-node2 Default=NO MaxTime=48:00:00 State=UP
PartitionName=normal Nodes=cc1,cc2,worker[1-3] Default=YES MaxTime=24:00:00 State=UP
PartitionName=bigmem Nodes=cr5,txt455 Default=NO MaxTime=72:00:00 State=UP
PartitionName=mixed Nodes=c1,cc1,cr5,txt455 Default=NO MaxTime=24:00:00 State=UP
PartitionName=all Nodes=c1,cc1,cc2,cr5,txt455,worker[1-3],gpu-node1,gpu-node2 Default=NO MaxTime=INFINITE State=UP
```

#### Node Naming Flexibility Options:

**1. Individual Node Definitions:**
```ini
# Each node defined separately
NodeName=server1 CPUs=4 RealMemory=8000 State=UNKNOWN
NodeName=workstation-lab2 CPUs=8 RealMemory=16000 State=UNKNOWN
NodeName=hpc-node-A CPUs=16 RealMemory=32000 State=UNKNOWN
```

**2. Comma-Separated Lists:**
```ini
# Group nodes with identical specifications
NodeName=node1,node2,server3,lab-pc5 CPUs=8 RealMemory=16000 State=UNKNOWN
NodeName=gpu1,gpu-server,ai-box CPUs=12 RealMemory=24000 Gres=gpu:1 State=UNKNOWN
```

**3. Mixed Patterns:**
```ini
# Combine ranges and individual names
NodeName=compute[1-5],special-node,backup-server CPUs=8 RealMemory=16000 State=UNKNOWN
```

**4. Complex Naming with Features:**
```ini
# Different node types with descriptive names
NodeName=intel-node1,intel-node2 CPUs=8 RealMemory=16000 Features=intel State=UNKNOWN
NodeName=amd-server1,amd-server2 CPUs=16 RealMemory=32000 Features=amd State=UNKNOWN
NodeName=gpu-workstation1 CPUs=12 RealMemory=24000 Gres=gpu:2 Features=nvidia State=UNKNOWN
NodeName=old-server,legacy-box CPUs=4 RealMemory=8000 Features=legacy State=UNKNOWN
```

#### Real-World Example Configuration:

```ini
# Control node
ControlMachine=slurm-master
ControlAddr=192.168.1.10

# Different types of compute nodes
NodeName=c1 NodeAddr=192.168.1.101 CPUs=8 RealMemory=16000 Gres=gpu:1 Features=gpu,tesla State=UNKNOWN
NodeName=cc1 NodeAddr=192.168.1.102 CPUs=4 RealMemory=8000 Features=basic State=UNKNOWN
NodeName=cc2 NodeAddr=192.168.1.103 CPUs=4 RealMemory=8000 Features=basic State=UNKNOWN
NodeName=cr5 NodeAddr=192.168.1.105 CPUs=16 RealMemory=64000 Features=bigmem State=UNKNOWN
NodeName=txt455 NodeAddr=192.168.1.201 CPUs=8 RealMemory=12000 Features=mixed State=UNKNOWN

# Additional nodes with different naming
NodeName=lab-pc1,lab-pc2 CPUs=6 RealMemory=12000 Features=lab State=UNKNOWN
NodeName=server-bio NodeAddr=192.168.1.150 CPUs=20 RealMemory=128000 Features=bigmem,bio State=UNKNOWN

# Partition definitions
PartitionName=gpu Nodes=c1 Default=NO MaxTime=48:00:00 State=UP AllowGroups=gpu-users
PartitionName=normal Nodes=cc1,cc2,txt455,lab-pc1,lab-pc2 Default=YES MaxTime=24:00:00 State=UP
PartitionName=bigmem Nodes=cr5,server-bio Default=NO MaxTime=72:00:00 State=UP AllowGroups=bigmem-users
PartitionName=bio Nodes=server-bio,txt455 Default=NO MaxTime=168:00:00 State=UP AllowGroups=biolab
PartitionName=all Nodes=c1,cc1,cc2,cr5,txt455,lab-pc1,lab-pc2,server-bio Default=NO MaxTime=INFINITE State=UP
```

#### Best Practices for Mixed Node Names:

**1. Use Descriptive Names:**
```ini
# Good: descriptive and meaningful
NodeName=gpu-workstation1,cpu-server1,bigmem-node1

# Avoid: cryptic names without context
NodeName=n1,x2,z99
```

**2. Include Node Addresses for Clarity:**
```ini
NodeName=c1 NodeAddr=192.168.1.101 CPUs=8 State=UNKNOWN
NodeName=cc1 NodeAddr=192.168.1.102 CPUs=4 State=UNKNOWN
```

**3. Use Features for Node Categorization:**
```ini
NodeName=c1 Features=gpu,nvidia CPUs=8 Gres=gpu:1 State=UNKNOWN
NodeName=cc1,cc2 Features=basic,intel CPUs=4 State=UNKNOWN
NodeName=cr5 Features=bigmem,amd CPUs=16 RealMemory=64000 State=UNKNOWN
```

**4. Group Similar Nodes When Possible:**
```ini
# Instead of repeating specifications
NodeName=cc1 CPUs=4 RealMemory=8000 State=UNKNOWN
NodeName=cc2 CPUs=4 RealMemory=8000 State=UNKNOWN

# Use grouped definition
NodeName=cc1,cc2 CPUs=4 RealMemory=8000 State=UNKNOWN
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

### Restarting SLURM Services

When you make configuration changes to SLURM, you'll need to restart the services. Here's how to do it properly:

#### Method 1: Restart Individual Services

**On the control node:**
```bash
# Restart the SLURM controller daemon
sudo systemctl restart slurmctld

# Check status
sudo systemctl status slurmctld
```

**On compute nodes:**
```bash
# Restart the SLURM daemon on each compute node
sudo systemctl restart slurmd

# Check status
sudo systemctl status slurmd
```

**On all nodes (MUNGE authentication):**
```bash
# If you changed MUNGE keys, restart MUNGE first
sudo systemctl restart munge
sudo systemctl status munge
```

#### Method 2: Graceful Restart with Configuration Reload

**Reload configuration without full restart:**
```bash
# On control node - reload configuration without stopping jobs
sudo scontrol reconfigure

# This tells all nodes to reload slurm.conf without restarting services
```

#### Method 3: Complete Cluster Restart

**For major configuration changes, restart in this order:**

1. **Stop all services (reverse order):**
```bash
# On compute nodes first
sudo systemctl stop slurmd

# On control node
sudo systemctl stop slurmctld

# On all nodes
sudo systemctl stop munge
```

2. **Start all services (correct order):**
```bash
# On all nodes - start MUNGE first
sudo systemctl start munge

# On control node - start controller
sudo systemctl start slurmctld

# On compute nodes - start worker daemons
sudo systemctl start slurmd
```

#### Method 4: Using SSH for Multi-Node Restart

**Restart all compute nodes from control node:**
```bash
# Create list of compute nodes
COMPUTE_NODES="c1 cc1 cc2 cr5 txt455"

# Restart slurmd on all compute nodes
for node in $COMPUTE_NODES; do
    echo "Restarting slurmd on $node"
    ssh $node "sudo systemctl restart slurmd"
done

# Check status on all nodes
for node in $COMPUTE_NODES; do
    echo "Status on $node:"
    ssh $node "sudo systemctl status slurmd --no-pager -l"
done
```

#### Method 5: Restart with Job Consideration

**Check for running jobs before restart:**
```bash
# Check if jobs are running
squeue

# Check node states
sinfo

# Drain nodes before restart (prevents new jobs from starting)
sudo scontrol update NodeName=c1,cc1,cc2 State=DRAIN Reason="Maintenance restart"

# Wait for jobs to finish or cancel them
# scancel <job_id>

# Restart services
sudo systemctl restart slurmd  # on compute nodes
sudo systemctl restart slurmctld  # on control node

# Resume nodes after restart
sudo scontrol update NodeName=c1,cc1,cc2 State=RESUME
```

#### Troubleshooting Restart Issues

**Common restart problems and solutions:**

1. **Check service status:**
```bash
sudo systemctl status slurmctld
sudo systemctl status slurmd
sudo systemctl status munge

# View recent logs
sudo journalctl -u slurmctld -f
sudo journalctl -u slurmd -f
```

2. **Check configuration syntax:**
```bash
# Test SLURM configuration
sudo slurmd -D  # Test on compute node
sudo slurmctld -D  # Test on control node
```

3. **Check MUNGE authentication:**
```bash
# Test MUNGE between nodes
munge -n | ssh compute-node unmunge
```

4. **Reset node states if needed:**
```bash
# If nodes show as DOWN after restart
sudo scontrol update NodeName=c1,cc1,cc2 State=IDLE
sudo scontrol update NodeName=ALL State=RESUME
```

### Troubleshooting SLURM Controller (slurmctld) Failures

When `slurmctld.service` fails to start, follow these diagnostic steps:

#### Step 1: Check Detailed Error Logs

```bash
# View detailed systemd logs
sudo journalctl -u slurmctld.service -n 50

# View SLURM-specific logs
sudo tail -f /var/log/slurm-llnl/slurmctld.log

# Check for any permission issues
ls -la /var/log/slurm-llnl/
```

#### Step 2: Test Configuration Manually

```bash
# Test slurmctld configuration in debug mode
sudo slurmctld -D

# This will show detailed error messages about configuration issues
# Press Ctrl+C to stop after seeing the errors
```

#### Step 3: Common Configuration Issues

**1. Check slurm.conf syntax:**
```bash
# Verify configuration file exists and has correct permissions
sudo ls -la /etc/slurm-llnl/slurm.conf

# Check for syntax errors
sudo slurmd -t  # Test configuration
```

**2. Verify MUNGE is running:**
```bash
# MUNGE must be running before slurmctld
sudo systemctl status munge
sudo systemctl start munge  # if not running
```

**3. Check hostname resolution:**
```bash
# Ensure ControlMachine hostname resolves correctly
hostname
hostname -f
ping $(hostname)

# Update /etc/hosts if needed
echo "127.0.0.1 $(hostname)" | sudo tee -a /etc/hosts
```

**4. Verify file permissions:**
```bash
# Check SLURM directories exist with correct permissions
sudo mkdir -p /var/log/slurm-llnl
sudo mkdir -p /var/spool/slurm
sudo chown slurm:slurm /var/log/slurm-llnl
sudo chown slurm:slurm /var/spool/slurm
```

#### Step 4: Fix Common Issues

**Issue 1: Configuration file not found**
```bash
# Error: "slurm.conf not found"
# Solution: Ensure file exists in correct location
sudo find / -name "slurm.conf" 2>/dev/null
sudo cp /path/to/slurm.conf /etc/slurm-llnl/slurm.conf
```

**Issue 2: MUNGE authentication failure**
```bash
# Error: "MUNGE authentication error"
# Solution: Restart MUNGE and check key
sudo systemctl restart munge
sudo systemctl status munge

# Test MUNGE locally
munge -n | unmunge
```

**Issue 3: Hostname mismatch**
```bash
# Error: "ControlMachine hostname mismatch"
# Solution: Update slurm.conf with correct hostname
sudo sed -i "s/ControlMachine=.*/ControlMachine=$(hostname)/" /etc/slurm-llnl/slurm.conf
```

**Issue 4: Port already in use**
```bash
# Error: "Port 6817 already in use"
# Solution: Check what's using the port
sudo netstat -tlnp | grep 6817
sudo lsof -i :6817

# Kill conflicting process or change port in slurm.conf
```

**Issue 5: Database connection failure (if accounting enabled)**
```bash
# Error: "Database connection failed"
# Solution: Check slurmdbd service and database
sudo systemctl status slurmdbd
sudo systemctl restart slurmdbd

# Test database connection
mysql -u slurm -p slurm_acct_db
```

#### Step 5: Reset and Clean Start

If issues persist, try a clean restart:

```bash
# Stop all SLURM services
sudo systemctl stop slurmctld
sudo systemctl stop slurmdbd  # if using accounting
sudo pkill -f slurm  # kill any remaining processes

# Clean temporary files
sudo rm -f /var/run/slurm*
sudo rm -f /tmp/slurm*

# Restart in correct order
sudo systemctl start munge
sudo systemctl start slurmdbd  # if using accounting
sudo systemctl start slurmctld

# Check status
sudo systemctl status slurmctld
```

#### Step 6: Minimal Configuration Test

Create a minimal working configuration for testing:

```bash
# Backup current config
sudo cp /etc/slurm-llnl/slurm.conf /etc/slurm-llnl/slurm.conf.backup

# Create minimal test config
sudo tee /etc/slurm-llnl/slurm.conf << EOF
ControlMachine=$(hostname)
ControlAddr=$(hostname -I | awk '{print $1}')
SlurmUser=slurm
SlurmdUser=root
StateSaveLocation=/var/spool/slurm
SlurmdSpoolDir=/var/spool/slurm/d
SwitchType=switch/none
MpiDefault=none
ProctrackType=proctrack/pgid
ReturnToService=1
SlurmctldPidFile=/var/run/slurmctld.pid
SlurmdPidFile=/var/run/slurmd.pid
SlurmdLogFile=/var/log/slurm-llnl/slurmd.log
SlurmctldLogFile=/var/log/slurm-llnl/slurmctld.log

# Node and partition (adjust as needed)
NodeName=$(hostname) CPUs=1 State=UNKNOWN
PartitionName=debug Nodes=$(hostname) Default=YES MaxTime=INFINITE State=UP
EOF

# Test with minimal config
sudo slurmctld -D
```

#### Step 7: Debugging Commands Summary

Quick diagnostic commands to run:

```bash
# Check all SLURM-related processes
ps aux | grep slurm

# Check SLURM configuration
sudo slurmd -t

# View configuration being used
sudo slurmctld -D 2>&1 | head -20

# Check network ports
sudo netstat -tlnp | grep -E "(6817|6818|6819)"

# Check SLURM user exists
id slurm

# Verify MUNGE key permissions
sudo ls -la /etc/munge/munge.key

# Check systemd service file
sudo systemctl cat slurmctld.service
```

#### Automated Restart Script

**Create a restart script for convenience:**
```bash
cat > /usr/local/bin/restart_slurm.sh << 'EOF'
#!/bin/bash

# SLURM restart script
COMPUTE_NODES="c1 cc1 cc2 cr5 txt455"

echo "Restarting SLURM cluster..."

# Function to restart on remote nodes
restart_compute_nodes() {
    for node in $COMPUTE_NODES; do
        echo "Restarting slurmd on $node"
        ssh $node "sudo systemctl restart munge && sudo systemctl restart slurmd"
        if [ $? -eq 0 ]; then
            echo "✓ $node restarted successfully"
        else
            echo "✗ Failed to restart $node"
        fi
    done
}

# Restart control node services
echo "Restarting control node services..."
sudo systemctl restart munge
sudo systemctl restart slurmctld

# Restart compute nodes
restart_compute_nodes

# Check cluster status
echo "Checking cluster status..."
sleep 5
sinfo
squeue

echo "SLURM restart completed!"
EOF

chmod +x /usr/local/bin/restart_slurm.sh
```

#### When to Restart vs Reconfigure

**Use `scontrol reconfigure` when:**
- Changing partition configurations
- Modifying node properties
- Updating time limits
- Adding/removing nodes

**Use full restart when:**
- Changing authentication settings
- Modifying core daemon configurations
- Updating SLURM software
- Changing database connections
- Installing new SLURM versions

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

SLURM accounting provides detailed tracking of job usage, resource consumption, and user activity. This is essential for resource allocation, billing, and cluster management.

### 6.1 Install and Configure MariaDB/MySQL

#### Option A: Traditional Installation

On the control node (or a dedicated database server):

```bash
# Install MariaDB
sudo apt update
sudo apt install mariadb-server mariadb-client
# Or for CentOS/RHEL:
# sudo yum install mariadb-server mariadb

# Secure the installation
sudo mysql_secure_installation
```

Create the SLURM accounting database:

```bash
sudo mysql -u root -p
```

```sql
CREATE DATABASE slurm_acct_db;
CREATE USER 'slurm'@'localhost' IDENTIFIED BY 'your_secure_password';
GRANT ALL PRIVILEGES ON slurm_acct_db.* TO 'slurm'@'localhost';
FLUSH PRIVILEGES;
EXIT;
```

#### Option B: Docker Container Database

Using Docker for the database server provides better isolation, easier management, and consistent deployment:

**Install Docker (if not already installed):**

```bash
# Install Docker
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh
sudo usermod -aG docker $USER
sudo systemctl enable docker
sudo systemctl start docker

# Install Docker Compose
sudo apt install docker-compose
```

**Create Docker Compose configuration:**

Create `/opt/slurm-db/docker-compose.yml`:

```yaml
version: '3.8'
services:
  mariadb:
    image: mariadb:10.11
    container_name: slurm-mariadb
    restart: unless-stopped
    environment:
      MYSQL_ROOT_PASSWORD: secure_root_password
      MYSQL_DATABASE: slurm_acct_db
      MYSQL_USER: slurm
      MYSQL_PASSWORD: your_secure_password
    ports:
      - "3306:3306"
    volumes:
      - mariadb_data:/var/lib/mysql
      - ./mariadb_config:/etc/mysql/conf.d
      - ./backups:/backups
    networks:
      - slurm-network

volumes:
  mariadb_data:

networks:
  slurm-network:
    driver: bridge
```

**Create MariaDB configuration file:**

Create `/opt/slurm-db/mariadb_config/slurm.cnf`:

```ini
[mysqld]
# SLURM specific optimizations
innodb_buffer_pool_size = 128M
innodb_log_file_size = 64M
max_connections = 100
query_cache_size = 32M
query_cache_type = 1

# Character set
character-set-server = utf8mb4
collation-server = utf8mb4_unicode_ci

# Binary logging (for backups)
log-bin = mysql-bin
binlog_format = ROW
expire_logs_days = 7
```

**Start the database container:**

```bash
# Create directories
sudo mkdir -p /opt/slurm-db/mariadb_config
sudo mkdir -p /opt/slurm-db/backups

# Create the configuration files (as shown above)

# Start the database
cd /opt/slurm-db
sudo docker-compose up -d

# Check if container is running
sudo docker-compose ps

# View logs
sudo docker-compose logs mariadb
```

**Connect to the containerized database:**

```bash
# Method 1: Using docker exec
sudo docker exec -it slurm-mariadb mysql -u root -p

# Method 2: Using mysql client on host (install mysql-client first)
sudo apt install mysql-client
mysql -h localhost -P 3306 -u slurm -p slurm_acct_db
```

**Additional Docker database management:**

```bash
# Stop the database
sudo docker-compose down

# Update the database (pull new image)
sudo docker-compose pull
sudo docker-compose up -d

# View container resource usage
sudo docker stats slurm-mariadb

# Backup database from container
sudo docker exec slurm-mariadb mysqldump -u slurm -p'your_secure_password' slurm_acct_db > /opt/slurm-db/backups/backup_$(date +%Y%m%d).sql

# Restore database to container
sudo docker exec -i slurm-mariadb mysql -u slurm -p'your_secure_password' slurm_acct_db < /opt/slurm-db/backups/backup_file.sql
```

### 6.2 Install and Configure SLURM Database Daemon (slurmdbd)

Install slurmdbd on the control node:

```bash
sudo apt install slurmdbd
# Or for CentOS/RHEL:
# sudo yum install slurm-slurmdbd
```

Create the slurmdbd configuration file `/etc/slurm-llnl/slurmdbd.conf`:

```ini
# slurmdbd.conf
AuthType=auth/munge
AuthInfo=/var/run/munge/munge.socket.2

# Database settings
DbdAddr=localhost
DbdHost=localhost
DbdPort=6819
SlurmUser=slurm

# Database connection
StorageType=accounting_storage/mysql
StorageHost=localhost
StoragePort=3306
StorageUser=slurm
StoragePass=your_secure_password
StorageLoc=slurm_acct_db

# Logging
LogFile=/var/log/slurm-llnl/slurmdbd.log
PidFile=/var/run/slurmdbd.pid
```

Set proper permissions:

```bash
sudo chown slurm:slurm /etc/slurm-llnl/slurmdbd.conf
sudo chmod 600 /etc/slurm-llnl/slurmdbd.conf
```

### 6.3 Update SLURM Configuration

**Note: This step is only for the control node (where slurmctld runs). Compute nodes do not need accounting configuration changes.**

Add accounting configuration to `/etc/slurm-llnl/slurm.conf` on the **control node only**:

```ini
# Add these lines to your existing slurm.conf on the CONTROL NODE
AccountingStorageType=accounting_storage/slurmdbd
AccountingStorageHost=localhost
AccountingStoragePort=6819
AccountingStorageEnforce=limits,safe
JobCompType=jobcomp/mysql
JobCompHost=localhost
JobCompUser=slurm
JobCompPass=your_secure_password
JobCompLoc=slurm_acct_db

# Optional: Enable job completion logging
JobCompType=jobcomp/filetxt
JobCompLoc=/var/log/slurm-llnl/job_completions.log
```

**Important Notes:**
- Only modify `slurm.conf` on the control node for accounting
- The standard `slurm.conf` (with node and partition definitions) should still be copied to all nodes
- Accounting services (`slurmdbd`) only run on the control node
- Compute nodes automatically report to the accounting system through the control node

### 6.4 Start Services and Initialize Database

Start the database and slurmdbd:

```bash
# Start MariaDB
sudo systemctl enable mariadb
sudo systemctl start mariadb

# Start slurmdbd
sudo systemctl enable slurmdbd
sudo systemctl start slurmdbd

# Check status
sudo systemctl status slurmdbd
```

Restart SLURM services:

```bash
sudo systemctl restart slurmctld
sudo systemctl restart slurmd  # on compute nodes
```

### 6.5 Set Up Accounting Hierarchy

Create clusters, accounts, and users:

```bash
# Add the cluster to accounting
sacctmgr add cluster cluster_name

# Add accounts (projects/groups)
sacctmgr add account biolab Description="Bioinformatics Lab"
sacctmgr add account physics Description="Physics Department"

# Add users to accounts
sacctmgr add user john Account=biolab
sacctmgr add user jane Account=physics

# Set resource limits (optional)
sacctmgr modify account biolab set MaxJobs=50 MaxNodes=10
sacctmgr modify user john set MaxJobs=20 MaxNodes=5
```

### 6.6 Useful Accounting Commands

Monitor and query accounting data:

```bash
# Show accounting configuration
sacctmgr show cluster
sacctmgr show account
sacctmgr show user

# Job accounting reports
sacct                           # Recent jobs
sacct -u username              # Jobs by user
sacct -A account_name          # Jobs by account
sacct -S 2025-01-01            # Jobs since date
sacct --format=JobID,JobName,Partition,Account,AllocCPUS,State,ExitCode

# Resource usage reports
sreport job sizes              # Job size distribution
sreport user top               # Top users by CPU hours
sreport cluster utilization    # Cluster utilization
sreport cluster AccountUtilizationByUser Start=2025-01-01
```

### 6.7 Fair Share Scheduling (Optional)

Enable fair share to prioritize jobs based on historical usage:

Add to `slurm.conf`:

```ini
PriorityType=priority/multifactor
PriorityDecayHalfLife=7-0      # 7 days
PriorityUsageResetPeriod=NONE
PriorityWeightFairshare=10000
PriorityWeightAge=1000
PriorityWeightPartition=1000
PriorityWeightJobSize=1000
PriorityMaxAge=7-0
```

### 6.8 Troubleshooting Accounting

Common issues and solutions:

```bash
# Check slurmdbd logs
sudo tail -f /var/log/slurm-llnl/slurmdbd.log

# Test database connection
mysql -u slurm -p slurm_acct_db

# Reset accounting if needed
sudo systemctl stop slurmdbd
mysql -u root -p -e "DROP DATABASE slurm_acct_db; CREATE DATABASE slurm_acct_db;"
sudo systemctl start slurmdbd

# Force accounting update
sudo scontrol reconfigure
```

### 6.9 Backup and Maintenance

#### Traditional Database Backup:

```bash
# Create backup script
cat > /usr/local/bin/backup_slurm_db.sh << 'EOF'
#!/bin/bash
BACKUP_DIR="/backup/slurm"
DATE=$(date +%Y%m%d_%H%M%S)
mkdir -p $BACKUP_DIR
mysqldump -u slurm -p'your_secure_password' slurm_acct_db > $BACKUP_DIR/slurm_acct_db_$DATE.sql
find $BACKUP_DIR -name "*.sql" -mtime +30 -delete
EOF

chmod +x /usr/local/bin/backup_slurm_db.sh

# Add to crontab for daily backups
echo "0 2 * * * /usr/local/bin/backup_slurm_db.sh" | sudo crontab -
```

#### Docker Database Backup:

```bash
# Create Docker-specific backup script
cat > /usr/local/bin/backup_slurm_docker_db.sh << 'EOF'
#!/bin/bash
BACKUP_DIR="/opt/slurm-db/backups"
DATE=$(date +%Y%m%d_%H%M%S)
CONTAINER_NAME="slurm-mariadb"

mkdir -p $BACKUP_DIR

# Method 1: Database dump from container
docker exec $CONTAINER_NAME mysqldump -u slurm -p'your_secure_password' slurm_acct_db > $BACKUP_DIR/slurm_acct_db_$DATE.sql

# Method 2: Full volume backup (includes all databases and configuration)
docker run --rm -v slurm-db_mariadb_data:/source -v $BACKUP_DIR:/backup alpine tar czf /backup/mariadb_volume_$DATE.tar.gz -C /source .

# Cleanup old backups (keep 30 days)
find $BACKUP_DIR -name "*.sql" -mtime +30 -delete
find $BACKUP_DIR -name "*.tar.gz" -mtime +30 -delete

# Log backup completion
echo "$(date): SLURM database backup completed" >> /var/log/slurm_backup.log
EOF

chmod +x /usr/local/bin/backup_slurm_docker_db.sh

# Add to crontab for daily backups
echo "0 2 * * * /usr/local/bin/backup_slurm_docker_db.sh" | sudo crontab -
```

#### Docker Database Monitoring:

```bash
# Create monitoring script for Docker database
cat > /usr/local/bin/monitor_slurm_docker_db.sh << 'EOF'
#!/bin/bash
CONTAINER_NAME="slurm-mariadb"

# Check if container is running
if ! docker ps | grep -q $CONTAINER_NAME; then
    echo "$(date): SLURM MariaDB container is not running!" >> /var/log/slurm_db_monitor.log
    # Optionally restart the container
    cd /opt/slurm-db && docker-compose up -d
fi

# Check database connectivity
if ! docker exec $CONTAINER_NAME mysqladmin -u slurm -p'your_secure_password' ping > /dev/null 2>&1; then
    echo "$(date): SLURM database is not responding!" >> /var/log/slurm_db_monitor.log
fi

# Check disk space
DISK_USAGE=$(docker exec $CONTAINER_NAME df -h /var/lib/mysql | awk 'NR==2 {print $5}' | sed 's/%//')
if [ "$DISK_USAGE" -gt 85 ]; then
    echo "$(date): SLURM database disk usage is $DISK_USAGE%" >> /var/log/slurm_db_monitor.log
fi
EOF

chmod +x /usr/local/bin/monitor_slurm_docker_db.sh

# Add to crontab for monitoring every 5 minutes
echo "*/5 * * * * /usr/local/bin/monitor_slurm_docker_db.sh" | sudo crontab -
```

#### Docker Database Maintenance:

```bash
# Update MariaDB container
cd /opt/slurm-db
sudo docker-compose pull mariadb
sudo docker-compose up -d

# View database logs
sudo docker-compose logs -f mariadb

# Database optimization (run monthly)
sudo docker exec slurm-mariadb mysql -u slurm -p'your_secure_password' slurm_acct_db -e "OPTIMIZE TABLE job_table, step_table, assoc_table;"

# Clean up old Docker images
sudo docker image prune -f
```

## References

- [SLURM Quick Start Guide](https://slurm.schedmd.com/quickstart.html)
- [SLURM Documentation](https://slurm.schedmd.com/documentation.html)
