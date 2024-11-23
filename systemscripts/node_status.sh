#!/bin/bash

# List of nodes
nodes=(node10 node11 node12 node13 node14 node15 node16 node20 node21 node22 node23 node24 node33 node34 node36 node37)

# Function to get node status
get_node_status() {
    node=$1
    ssh $node "echo -n '$node: '; uptime | awk '{print \$10, \$11, \$12}'; echo -n 'Memory: '; free -h | awk '/^Mem:/ {print \$3\"/\"\$2}'; echo -n 'CPU: '; top -bn1 | grep '%Cpu' | awk '{print \$2\"%\"}'" 2>/dev/null
}

# Print header
printf "%-10s %-20s %-15s %-10s\n" "Node" "Load Average" "Memory" "CPU Usage"
echo "------------------------------------------------------"

# Loop through nodes and get status
for node in "${nodes[@]}"; do
    status=$(get_node_status $node)
    if [ ! -z "$status" ]; then
        echo "$status"
    else
        echo "$node: Unavailable"
    fi
done
