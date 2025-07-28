# Building Multi-Platform Docker Images

This guide explains how to build Docker images that support multiple CPU architectures (AMD64, ARM64) using GitHub Actions. Multi-platform images allow your containers to run on different hardware architectures seamlessly.

## Overview

Multi-platform Docker images enable your application to run on:
- **AMD64 (x86_64)**: Traditional Intel/AMD processors (most servers, laptops)
- **ARM64 (aarch64)**: ARM processors (Apple Silicon, AWS Graviton, Raspberry Pi)
- **Other architectures**: ARM v7, s390x, ppc64le (as needed)

## Architecture Strategies

### Strategy 1: Native Builds (Recommended for CI/CD)

Build each architecture on native hardware for optimal performance and compatibility:

- **AMD64**: Use standard `ubuntu-latest` runners
- **ARM64**: Use `ubuntu-22.04-arm` or ARM-based runners
- **Merge**: Combine into a single multi-platform manifest

### Strategy 2: Cross-Compilation with QEMU

Build all architectures on a single machine using emulation:

- **Pros**: Simpler setup, single runner
- **Cons**: Slower builds, potential compatibility issues

## Method 1: Native Multi-Architecture Builds

This approach builds each architecture on native hardware for best performance.

### Prerequisites

1. **GitHub Secrets Setup**:
   ```
   DH_USER: Your Docker Hub username
   DH_TOKEN: Your Docker Hub access token
   ```

2. **Runner Requirements**:
   - Standard GitHub runners for AMD64
   - ARM64 runners (GitHub-hosted or self-hosted)

### GitHub Actions Workflow

#### Step 1: AMD64 Build Job

```yaml
name: Multi-Platform Docker Build

on:
  push:
    branches: [ main ]
    paths:
      - 'your-app/**'
  release:
    types: [published]
  workflow_dispatch:

env:
  REGISTRY_IMAGE: your-dockerhub-user/your-app

jobs:
  amd64:
    runs-on: ubuntu-latest
    steps:
      - name: Checkout repository
        uses: actions/checkout@v4

      - name: Set up Docker Buildx
        uses: docker/setup-buildx-action@v3

      - name: Login to Docker Hub
        uses: docker/login-action@v3
        with:
          username: ${{ secrets.DH_USER }}
          password: ${{ secrets.DH_TOKEN }}

      - name: Extract metadata
        id: meta
        uses: docker/metadata-action@v5
        with:
          images: ${{ env.REGISTRY_IMAGE }}

      - name: Build and push AMD64
        id: build
        uses: docker/build-push-action@v5
        with:
          context: ./your-app/
          platforms: linux/amd64
          labels: ${{ steps.meta.outputs.labels }}
          cache-from: type=registry,ref=${{ env.REGISTRY_IMAGE }}:amd64-cache
          cache-to: type=registry,ref=${{ env.REGISTRY_IMAGE }}:amd64-cache,mode=max
          outputs: type=image,push-by-digest=true,name-canonical=true,push=true

      - name: Export digest
        run: |
          mkdir -p ${{ runner.temp }}/digests
          digest="${{ steps.build.outputs.digest }}"
          touch "${{ runner.temp }}/digests/${digest#sha256:}"

      - name: Upload digest
        uses: actions/upload-artifact@v4
        with:
          name: digests-amd64
          path: ${{ runner.temp }}/digests/*
          if-no-files-found: error
          retention-days: 1
```

#### Step 2: ARM64 Build Job

```yaml
  arm64:
    runs-on: ubuntu-22.04-arm  # or your ARM64 runner
    steps:
      - name: Checkout repository
        uses: actions/checkout@v4

      - name: Set up Docker Buildx
        uses: docker/setup-buildx-action@v3

      - name: Login to Docker Hub
        uses: docker/login-action@v3
        with:
          username: ${{ secrets.DH_USER }}
          password: ${{ secrets.DH_TOKEN }}

      - name: Extract metadata
        id: meta
        uses: docker/metadata-action@v5
        with:
          images: ${{ env.REGISTRY_IMAGE }}

      - name: Build and push ARM64
        id: build
        uses: docker/build-push-action@v5
        with:
          context: ./your-app/
          platforms: linux/arm64
          labels: ${{ steps.meta.outputs.labels }}
          cache-from: type=registry,ref=${{ env.REGISTRY_IMAGE }}:arm64-cache
          cache-to: type=registry,ref=${{ env.REGISTRY_IMAGE }}:arm64-cache,mode=max
          outputs: type=image,push-by-digest=true,name-canonical=true,push=true

      - name: Export digest
        run: |
          mkdir -p ${{ runner.temp }}/digests
          digest="${{ steps.build.outputs.digest }}"
          touch "${{ runner.temp }}/digests/${digest#sha256:}"

      - name: Upload digest
        uses: actions/upload-artifact@v4
        with:
          name: digests-arm64
          path: ${{ runner.temp }}/digests/*
          if-no-files-found: error
          retention-days: 1
```

#### Step 3: Merge and Create Manifest

```yaml
  merge:
    runs-on: ubuntu-latest
    needs:
      - amd64
      - arm64
    steps:
      - name: Download digests
        uses: actions/download-artifact@v4
        with:
          path: ${{ runner.temp }}/digests
          pattern: digests-*
          merge-multiple: true

      - name: Set up Docker Buildx
        uses: docker/setup-buildx-action@v3

      - name: Login to Docker Hub
        uses: docker/login-action@v3
        with:
          username: ${{ secrets.DH_USER }}
          password: ${{ secrets.DH_TOKEN }}

      - name: Extract metadata
        id: meta
        uses: docker/metadata-action@v5
        with:
          images: ${{ env.REGISTRY_IMAGE }}
          tags: |
            type=ref,event=branch
            type=ref,event=pr
            type=semver,pattern={{version}}
            type=semver,pattern={{major}}.{{minor}}
            type=raw,value=latest,enable={{is_default_branch}}

      - name: Create manifest list and push
        working-directory: ${{ runner.temp }}/digests
        run: |
          docker buildx imagetools create $(jq -cr '.tags | map("-t " + .) | join(" ")' <<< "$DOCKER_METADATA_OUTPUT_JSON") \
            $(printf '${{ env.REGISTRY_IMAGE }}@sha256:%s ' *)

      - name: Inspect image
        run: |
          docker buildx imagetools inspect ${{ env.REGISTRY_IMAGE }}:${{ steps.meta.outputs.version }}
```

## Method 2: QEMU Cross-Compilation (Single Runner)

For simpler setups, you can build all architectures on a single runner using QEMU emulation:

```yaml
name: Multi-Platform Docker Build (QEMU)

on:
  push:
    branches: [ main ]
  workflow_dispatch:

env:
  REGISTRY_IMAGE: your-dockerhub-user/your-app

jobs:
  build:
    runs-on: ubuntu-latest
    steps:
      - name: Checkout
        uses: actions/checkout@v4

      - name: Set up QEMU
        uses: docker/setup-qemu-action@v3

      - name: Set up Docker Buildx
        uses: docker/setup-buildx-action@v3

      - name: Login to Docker Hub
        uses: docker/login-action@v3
        with:
          username: ${{ secrets.DH_USER }}
          password: ${{ secrets.DH_TOKEN }}

      - name: Extract metadata
        id: meta
        uses: docker/metadata-action@v5
        with:
          images: ${{ env.REGISTRY_IMAGE }}
          tags: |
            type=ref,event=branch
            type=semver,pattern={{version}}
            type=raw,value=latest,enable={{is_default_branch}}

      - name: Build and push multi-platform
        uses: docker/build-push-action@v5
        with:
          context: .
          platforms: linux/amd64,linux/arm64
          push: true
          tags: ${{ steps.meta.outputs.tags }}
          labels: ${{ steps.meta.outputs.labels }}
          cache-from: type=gha
          cache-to: type=gha,mode=max
```

## Docker Configuration Optimization

### Docker Daemon Configuration

Create `daemon.json` for optimized builds:

```json
{
  "experimental": true,
  "features": {
    "buildkit": true
  },
  "data-root": "/mnt/docker",
  "storage-driver": "overlay2",
  "storage-opts": [
    "overlay2.override_kernel_check=true"
  ],
  "registry-mirrors": [],
  "insecure-registries": [],
  "debug": false,
  "hosts": ["fd://"],
  "log-level": "info",
  "max-concurrent-downloads": 10,
  "max-concurrent-uploads": 5
}
```

### Apply Docker Configuration in CI:

```yaml
      - name: Configure Docker daemon
        run: |
          sudo mkdir -p /mnt/docker
          sudo chmod 777 /mnt/docker
          sudo mkdir -p /etc/docker
          sudo cp daemon.json /etc/docker/daemon.json
          sudo systemctl restart docker
          sudo docker info
```

## Advanced Configurations

### Custom Build Arguments per Architecture

```yaml
      - name: Build with architecture-specific args
        uses: docker/build-push-action@v5
        with:
          context: .
          platforms: linux/amd64,linux/arm64
          build-args: |
            ARCH=${{ matrix.arch }}
            VARIANT=${{ matrix.variant }}
          push: true
          tags: ${{ steps.meta.outputs.tags }}
```

### Matrix Strategy for Multiple Architectures

```yaml
jobs:
  build:
    runs-on: ${{ matrix.runner }}
    strategy:
      fail-fast: false
      matrix:
        include:
          - platform: linux/amd64
            runner: ubuntu-latest
            arch: amd64
          - platform: linux/arm64
            runner: ubuntu-22.04-arm
            arch: arm64
          - platform: linux/arm/v7
            runner: ubuntu-latest
            arch: armv7
    steps:
      - name: Build for ${{ matrix.platform }}
        uses: docker/build-push-action@v5
        with:
          platforms: ${{ matrix.platform }}
          # ... other configuration
```

### BuildKit Features

Enable advanced BuildKit features:

```yaml
      - name: Build with BuildKit features
        uses: docker/build-push-action@v5
        with:
          context: .
          platforms: linux/amd64,linux/arm64
          provenance: true
          sbom: true
          push: true
          tags: ${{ steps.meta.outputs.tags }}
          cache-from: |
            type=gha
            type=registry,ref=${{ env.REGISTRY_IMAGE }}:cache
          cache-to: |
            type=gha,mode=max
            type=registry,ref=${{ env.REGISTRY_IMAGE }}:cache,mode=max
```

## Dockerfile Best Practices for Multi-Platform

### Architecture-Aware Base Images

```dockerfile
# Use multi-platform base images
FROM node:18-alpine

# Or use architecture-specific logic
ARG TARGETARCH
RUN if [ "$TARGETARCH" = "arm64" ]; then \
      echo "ARM64 specific setup"; \
    elif [ "$TARGETARCH" = "amd64" ]; then \
      echo "AMD64 specific setup"; \
    fi
```

### Conditional Dependencies

```dockerfile
# Install architecture-specific packages
ARG TARGETARCH
RUN case "${TARGETARCH}" in \
      amd64) apt-get install -y x86-64-specific-package ;; \
      arm64) apt-get install -y arm64-specific-package ;; \
      *) echo "Unsupported architecture: ${TARGETARCH}" && exit 1 ;; \
    esac
```

### Multi-Stage Builds

```dockerfile
# Build stage
FROM --platform=$BUILDPLATFORM node:18-alpine AS builder
ARG TARGETARCH
ARG BUILDPLATFORM
WORKDIR /app
COPY package*.json ./
RUN npm ci --only=production

# Runtime stage
FROM node:18-alpine AS runtime
ARG TARGETARCH
COPY --from=builder /app/node_modules ./node_modules
COPY . .
EXPOSE 3000
CMD ["node", "server.js"]
```

## Testing Multi-Platform Images

### Local Testing

```bash
# Pull and test multi-platform image
docker pull your-dockerhub-user/your-app:latest

# Test on current architecture
docker run --rm your-dockerhub-user/your-app:latest

# Inspect manifest
docker buildx imagetools inspect your-dockerhub-user/your-app:latest
```

### CI/CD Testing

```yaml
  test:
    runs-on: ${{ matrix.runner }}
    needs: merge
    strategy:
      matrix:
        include:
          - runner: ubuntu-latest
            arch: amd64
          - runner: ubuntu-22.04-arm
            arch: arm64
    steps:
      - name: Test on ${{ matrix.arch }}
        run: |
          docker pull ${{ env.REGISTRY_IMAGE }}:latest
          docker run --rm ${{ env.REGISTRY_IMAGE }}:latest --version
```

## Troubleshooting

### Common Issues and Solutions

**1. QEMU Emulation Failures**
```bash
# Increase emulation timeout
docker buildx create --use --config buildkitd.toml
# buildkitd.toml:
# [worker.oci]
#   max-parallelism = 1
```

**2. Architecture-Specific Dependencies**
```dockerfile
# Use package manager conditionals
RUN apt-get update && \
    case "${TARGETARCH}" in \
      amd64) apt-get install -y intel-specific-package ;; \
      arm64) apt-get install -y arm-specific-package ;; \
    esac
```

**3. Build Cache Issues**
```yaml
# Use separate caches per architecture
cache-from: type=registry,ref=${{ env.REGISTRY_IMAGE }}:${{ matrix.arch }}-cache
cache-to: type=registry,ref=${{ env.REGISTRY_IMAGE }}:${{ matrix.arch }}-cache,mode=max
```

### Debugging Commands

```bash
# Check available platforms
docker buildx ls

# Inspect image manifest
docker manifest inspect your-image:tag

# Check platform-specific details
docker buildx imagetools inspect your-image:tag --format "{{json .}}"
```

## Performance Optimization

### Build Parallelization

```yaml
      - name: Build with parallel jobs
        uses: docker/build-push-action@v5
        with:
          context: .
          platforms: linux/amd64,linux/arm64
          build-args: |
            BUILDKIT_INLINE_CACHE=1
            MAKEFLAGS=-j$(nproc)
          cache-from: type=gha
          cache-to: type=gha,mode=max
```

### Registry Optimization

```yaml
      - name: Optimize registry pushes
        uses: docker/build-push-action@v5
        with:
          context: .
          platforms: linux/amd64,linux/arm64
          push: true
          tags: ${{ steps.meta.outputs.tags }}
          # Push by digest for better deduplication
          outputs: type=image,push-by-digest=true,name-canonical=true,push=true
```

## Security Considerations

### Image Scanning

```yaml
  security:
    runs-on: ubuntu-latest
    needs: merge
    steps:
      - name: Run Trivy vulnerability scanner
        uses: aquasecurity/trivy-action@master
        with:
          image-ref: ${{ env.REGISTRY_IMAGE }}:latest
          format: 'sarif'
          output: 'trivy-results.sarif'

      - name: Upload Trivy scan results
        uses: github/codeql-action/upload-sarif@v2
        if: always()
        with:
          sarif_file: 'trivy-results.sarif'
```

### SBOM Generation

```yaml
      - name: Build with SBOM
        uses: docker/build-push-action@v5
        with:
          context: .
          platforms: linux/amd64,linux/arm64
          push: true
          tags: ${{ steps.meta.outputs.tags }}
          sbom: true
          provenance: true
```

This comprehensive guide covers all aspects of building multi-platform Docker images, from basic setup to advanced optimization and security considerations.