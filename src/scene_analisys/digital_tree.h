

#ifndef DIGITAL_TREE_H_
#define DIGITAL_TREE_H_

namespace digital_tree {
	template<typename T>
	class Node {
		public:
			T* value;
			std::vector<Node<T>*> children;
			Node() {

			}

			Node(int chSize) {
				children.reserve(chSize);
			}

			Node(T* value, int chSize) {
				this->value = value;
				children.reserve(chSize);
			}

			virtual ~Node() {
				delete value;
				/* the children of the node are supposed to have been already destroyed */
			}

			void setNrChildren(int size) {
				children.resize(size);
			}
	};

	template<typename T>
	class Tree {
		private:
			Node<T> root;
			int maxIdx; /* the maximum index of a number that might appear in a configuration */

			void postorderDestr(Node<T>* crt) {
				for (int i = 0; i < crt->children.size(); i++) {
					if (crt->children[i] != NULL) {
						postorderDestr(crt->children[i]);
					}
				}
				delete crt;
			}

		public:

			Tree(int maxIdx) {
				this->maxIdx = maxIdx;
				root.setNrChildren(maxIdx);
			}

			virtual ~Tree() {
				/* do a postorder transversal on the tree and call the destuctor on the nodes */
				for (int i = 0; i < root.children.size(); i++) {
					postorderDestr(root.children[i]);
				}
			}

			/** remove all the nodes in the tree */
			void clear() {
				/* do a postorder transversal on the tree and call the destuctor on the nodes */
				for (int i = 0; i < root.children.size(); i++) {
					postorderDestr(root.children[i]);
				}
				root.children.clear();
			}

			void resize(int maxIdx) {
				this->maxIdx = maxIdx;
				root.setNrChildren(maxIdx);
			}

			/** returns the node stored on the given path */
			T* findNode(const std::vector<int>& config) {
				Node<T>* crt = &root;
				for (int i = 0; i < config.size(); i++) {
					if (crt->children[config[i]] == NULL) {
						return NULL;
					}
					crt = crt->children[config[i]];
				}
				return crt->value;
			}

			/** insert a node in the tree according to a given configuration */
			void insertNode(const std::vector<int>& config, const T* val) {
				int i;
				Node<T>* crt = &root;
				for (i = 0; i < config.size(); i++) {
					if (crt->children[config[i]] == NULL) {
						break;
					}
					crt = crt->children[config[i]];
				}
				if (i == config.size()) {
					/* update the existing node */
					*(crt->value) = *val;
					return;
				}
				for (; i < config.size(); i++) {
					crt->children[config[i]] = new Node<T> (maxIdx);
					crt = crt->children[config[i]];
				}
				crt->value = val;
			}
	};
}
#endif /* DIGITAL_TREE_H_ */
