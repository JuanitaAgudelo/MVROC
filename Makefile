.PHONY: help stash

help:
	@echo "MVROC"
	@echo ""
	@echo "  make stash  - Stash local changes, pull from remote, then restore the stash"
	@echo ""

# Same workflow as in seap-udea / DRZ site Makefiles: if the tree is dirty,
# stash → pull --no-rebase → stash pop; if clean, just pull.
stash:
	@echo "▶  Syncing with remote (stash → pull → stash pop)…"
	@if git diff-index --quiet HEAD --; then \
		git pull --no-rebase; \
	else \
		git stash && git pull --no-rebase && git stash pop || \
			echo "⚠️  Check for conflicts after restoring local changes."; \
	fi
