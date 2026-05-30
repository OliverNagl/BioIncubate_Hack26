# ── Sync data/ with the digs cluster ──────────────────────────────────────
# Remote: olnagl@digs:~/pipelines/t3_capsids/data  ("digs" = your ssh config host)
#
#   make push        # local data/  ->  cluster  (upload)
#   make pull        # cluster      ->  local data/  (download)
#   make push-dry    # preview upload (rsync --dry-run, no changes)
#   make pull-dry    # preview download
#   make push-mirror # upload AND delete remote files missing locally (careful)
#   make pull-mirror # download AND delete local files missing remotely (careful)
#
# Override any variable on the command line, e.g.:
#   make pull REMOTE=olnagl@digs REMOTE_DIR='~/pipelines/t3_capsids/data'

REMOTE      ?= olnagl@digs
REMOTE_DIR  ?= ~/pipelines/t3_capsids/data
LOCAL_DIR   ?= data

# -a archive, -v verbose, -z compress, -h human sizes, --progress per-file bar
RSYNC       ?= rsync -avzh --progress

# Trailing slashes are deliberate: sync the *contents* of one dir into the other.
LOCAL       := $(LOCAL_DIR)/
REMOTE_PATH := $(REMOTE):$(REMOTE_DIR)/

.PHONY: help push pull push-dry pull-dry push-mirror pull-mirror remote-dir

help:
	@echo "Targets: push | pull | push-dry | pull-dry | push-mirror | pull-mirror"
	@echo "  REMOTE=$(REMOTE)  REMOTE_DIR=$(REMOTE_DIR)  LOCAL_DIR=$(LOCAL_DIR)"

# Ensure the remote directory exists before the first upload.
remote-dir:
	ssh $(REMOTE) 'mkdir -p $(REMOTE_DIR)'

push: remote-dir
	$(RSYNC) $(LOCAL) $(REMOTE_PATH)

pull:
	mkdir -p $(LOCAL_DIR)
	$(RSYNC) $(REMOTE_PATH) $(LOCAL)

push-dry: remote-dir
	$(RSYNC) --dry-run $(LOCAL) $(REMOTE_PATH)

pull-dry:
	$(RSYNC) --dry-run $(REMOTE_PATH) $(LOCAL)

push-mirror: remote-dir
	$(RSYNC) --delete $(LOCAL) $(REMOTE_PATH)

pull-mirror:
	mkdir -p $(LOCAL_DIR)
	$(RSYNC) --delete $(REMOTE_PATH) $(LOCAL)
